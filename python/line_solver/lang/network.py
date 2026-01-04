"""
Network class for LINE native Python implementation.

This module provides the Network class for constructing and managing queueing networks.
Ported from MATLAB implementation in matlab/src/lang/@MNetwork/
"""

from typing import Optional, Dict, List, Tuple, Union
import numpy as np

from .base import (
    Element, ElementType, Network as NetworkBase, NetworkElement, Node, StatefulNode, Station,
    JobClass, JobClassType, NodeType, SchedStrategy, RoutingStrategy
)
from .region import Region
from .routing import RoutingMatrix
from ..api.sn.network_struct import NetworkStruct
from ..constants import ProcessType


class Network(NetworkBase, Element):
    """
    Queueing network model for LINE.

    The Network class represents a complete queueing network model with nodes,
    job classes, and routing. It compiles to a NetworkStruct for solver consumption.
    """

    def __init__(self, name: str = "LineNetwork"):
        """
        Initialize a network.

        Args:
            name: Network name (default: "LineNetwork")
        """
        Element.__init__(self, ElementType.MODEL, name)
        self._nodes = []  # List of all nodes
        self._stations = []  # List of station nodes only
        self._classes = []  # List of job classes
        self._regions = []  # List of finite capacity regions
        self._connections = None  # Adjacency matrix (lazy-computed)
        self._routing_matrix = None  # Routing matrix
        self._sn = None  # Compiled NetworkStruct
        self._has_struct = False  # Whether struct is valid
        self._rewards = {}  # Dict[reward_name, reward_function]
        self._source_idx = -1  # Cached source node index
        self._sink_idx = -1  # Cached sink node index
        self.allow_replace = False  # When True, add_node replaces nodes with same name

    # =====================================================================
    # CORE CONSTRUCTION METHODS
    # =====================================================================

    def add_node(self, node: Node) -> None:
        """
        Add a node to the network.

        Args:
            node: Node to add (Queue, Source, Sink, Delay, Fork, Join, etc.)

        Raises:
            ValueError: If node already in network (when allow_replace is False)
        """
        if node in self._nodes:
            raise ValueError(f"[{self.name}] Node '{node.name}' already in network")

        # Check if replacing an existing node with the same name
        if self.allow_replace:
            for i, existing in enumerate(self._nodes):
                if existing.name == node.name:
                    # Replace existing node
                    node.set_model(self)
                    node._set_index(i)
                    self._nodes[i] = node

                    # Handle station replacement
                    if node.is_station():
                        old_station_idx = existing.get_station_index0()
                        if old_station_idx is not None and 0 <= old_station_idx < len(self._stations):
                            node._set_station_index(old_station_idx)
                            self._stations[old_station_idx] = node
                        else:
                            # Old node wasn't a station, add new one
                            node._set_station_index(len(self._stations))
                            self._stations.append(node)

                    self._reset_struct()
                    return

        # Normal case: add new node
        node.set_model(self)
        node._set_index(len(self._nodes))

        self._nodes.append(node)

        # Track stations separately
        if node.is_station():
            node._set_station_index(len(self._stations))
            self._stations.append(node)

        self._reset_struct()

    def add_class(self, jobclass: JobClass) -> None:
        """
        Add a job class to the network.

        Args:
            jobclass: Job class to add (OpenClass or ClosedClass)

        Raises:
            ValueError: If class already in network
        """
        if jobclass in self._classes:
            raise ValueError(f"[{self.name}] Class '{jobclass.name}' already in network")

        jobclass._set_index(len(self._classes))
        self._classes.append(jobclass)

        self._reset_struct()

    # MATLAB compatibility alias
    def addClass(self, jobclass: JobClass) -> None:
        """MATLAB-compatible alias for add_class()."""
        self.add_class(jobclass)

    def add_link(self, source: Node, dest: Node) -> None:
        """
        Add a link between two nodes.

        Args:
            source: Source node
            dest: Destination node

        Raises:
            ValueError: If nodes not in network
        """
        if source not in self._nodes:
            raise ValueError(f"[{self.name}] Source node '{source.name}' not in network")
        if dest not in self._nodes:
            raise ValueError(f"[{self.name}] Dest node '{dest.name}' not in network")

        # Initialize routing matrix if needed
        if self._routing_matrix is None:
            self._routing_matrix = RoutingMatrix(self)

        # Add route from source to dest for all classes with probability 1.0
        for job_class in self._classes:
            self._routing_matrix.set(job_class, job_class, source, dest, 1.0)

        # Invalidate connections and struct
        self._connections = None
        self._reset_struct()

    def add_links(self, routing_matrix: Union[RoutingMatrix, np.ndarray]) -> None:
        """
        Add multiple links via routing matrix.

        Args:
            routing_matrix: RoutingMatrix or numpy array with probabilities
        """
        if isinstance(routing_matrix, np.ndarray):
            routing_matrix = RoutingMatrix(routing_matrix)

        self._routing_matrix = routing_matrix
        self._connections = None
        self._reset_struct()

    def add_link_list(self, *nodes: Node) -> None:
        """
        Create serial links between nodes.

        Args:
            nodes: Variable number of nodes to link in sequence
        """
        for i in range(len(nodes) - 1):
            self.add_link(nodes[i], nodes[i + 1])

    def add_region(self, name_or_node: Union[str, Node], *nodes: Node) -> Region:
        """
        Add a finite capacity region containing the specified nodes.

        A finite capacity region constrains the total number of jobs
        that can be present in a group of nodes simultaneously.

        Args:
            name_or_node: Either a name for the region (str) or the first node.
                         If a string is provided, it's used as the region name.
                         If a node is provided, an auto-generated name is used.
            *nodes: Additional nodes to include in the region

        Returns:
            The created Region object for further configuration

        Examples:
            # With auto-generated name (FCR1, FCR2, etc.)
            fcr = model.add_region(queue1)
            fcr = model.add_region(queue1, queue2)

            # With custom name
            fcr = model.add_region('MyRegion', queue1)
            fcr = model.add_region('MyRegion', queue1, queue2)
        """
        if isinstance(name_or_node, str):
            # First argument is a name
            region_name = name_or_node
            node_list = list(nodes)
        else:
            # First argument is a node, auto-generate name
            region_name = f"FCR{len(self._regions) + 1}"
            node_list = [name_or_node] + list(nodes)

        region = Region(node_list, self._classes)
        region.set_name(region_name)
        self._regions.append(region)
        self._reset_struct()
        return region

    # MATLAB-style alias
    def addRegion(self, name_or_node: Union[str, Node], *nodes: Node) -> Region:
        """MATLAB-compatible alias for add_region()."""
        return self.add_region(name_or_node, *nodes)

    @property
    def regions(self) -> List[Region]:
        """Get the list of finite capacity regions."""
        return self._regions

    def get_regions(self) -> List[Region]:
        """Get the list of finite capacity regions."""
        return self._regions

    # MATLAB-style alias
    def getRegions(self) -> List[Region]:
        """MATLAB-compatible alias for get_regions()."""
        return self.get_regions()

    # =====================================================================
    # ROUTING METHODS
    # =====================================================================

    def init_routing_matrix(self) -> RoutingMatrix:
        """
        Initialize an empty routing matrix for all nodes and classes.

        Returns:
            RoutingMatrix with zeros (no routing initially)
        """
        return RoutingMatrix(self)

    def link(self, routing_matrix) -> None:
        """
        Set the routing matrix for the network.

        Args:
            routing_matrix: RoutingMatrix, dict with (from_class, to_class) keys,
                           or 2D list/array (applies same routing to all classes)

        Raises:
            ValueError: If routing matrix dimensions don't match network
        """
        if isinstance(routing_matrix, RoutingMatrix):
            self._routing_matrix = routing_matrix
        elif isinstance(routing_matrix, dict):
            # Convert dict to RoutingMatrix format
            # The dict has (from_class, to_class) tuple keys mapping to numpy arrays
            rm = self.init_routing_matrix()
            for (from_class, to_class), matrix in routing_matrix.items():
                matrix_arr = np.asarray(matrix)
                # Set routing probabilities for each node pair
                nodes = self.get_nodes()
                for i in range(min(len(nodes), matrix_arr.shape[0])):
                    for j in range(min(len(nodes), matrix_arr.shape[1])):
                        if matrix_arr[i, j] > 0:
                            rm.set(from_class, to_class, nodes[i], nodes[j], matrix_arr[i, j])
            self._routing_matrix = rm
        elif isinstance(routing_matrix, list):
            # Handle nested list format P[r][s] or simple 2D list
            rm = self.init_routing_matrix()
            nodes = self.get_nodes()
            classes = self.get_classes()
            nclasses = len(classes)

            # Check if it's a nested list P[r][s] format where each P[r][s] is a node routing matrix
            # P[r][s] format: routing_matrix[r][s] is a 2D array of node routing probabilities
            # This requires routing_matrix[r][s] to be a list/array, not a scalar
            is_class_routing = False
            if (len(routing_matrix) == nclasses and
                isinstance(routing_matrix[0], list) and
                len(routing_matrix[0]) == nclasses):
                # Check if routing_matrix[0][0] is a list/array (node routing matrix)
                # If it's a number, this is simple node routing, not class routing
                if isinstance(routing_matrix[0][0], (list, np.ndarray)):
                    is_class_routing = True

            if is_class_routing:
                # P[r][s] format - routing from class r to class s
                for r in range(nclasses):
                    for s in range(nclasses):
                        matrix_arr = np.asarray(routing_matrix[r][s])
                        if matrix_arr.size == 0:
                            continue
                        for i in range(min(len(nodes), matrix_arr.shape[0])):
                            for j in range(min(len(nodes), matrix_arr.shape[1])):
                                if matrix_arr[i, j] > 0:
                                    rm.set(classes[r], classes[s], nodes[i], nodes[j], matrix_arr[i, j])
            else:
                # Simple 2D list - apply same routing to all classes
                matrix_arr = np.asarray(routing_matrix)
                for jobclass in classes:
                    for i in range(min(len(nodes), matrix_arr.shape[0])):
                        for j in range(min(len(nodes), matrix_arr.shape[1])):
                            if matrix_arr[i, j] > 0:
                                rm.set(jobclass, jobclass, nodes[i], nodes[j], matrix_arr[i, j])
            self._routing_matrix = rm
        elif isinstance(routing_matrix, np.ndarray):
            # Convert numpy array to RoutingMatrix
            rm = self.init_routing_matrix()
            nodes = self.get_nodes()
            classes = self.get_classes()
            matrix_arr = routing_matrix
            for jobclass in classes:
                for i in range(min(len(nodes), matrix_arr.shape[0])):
                    for j in range(min(len(nodes), matrix_arr.shape[1])):
                        if matrix_arr[i, j] > 0:
                            rm.set(jobclass, jobclass, nodes[i], nodes[j], matrix_arr[i, j])
            self._routing_matrix = rm
        else:
            raise ValueError("[{0}] routing_matrix must be RoutingMatrix, dict, or 2D list/array".format(self.name))

        # Insert auto-generated ClassSwitch nodes for class switching routes
        self._insert_auto_class_switches()

        self._refresh_routing_matrix()
        self._reset_struct()

    def _insert_auto_class_switches(self) -> None:
        """
        Insert auto-generated ClassSwitch nodes for class switching in routing.

        When routing has class switching (jobs change class when moving between nodes),
        this method creates implicit ClassSwitch nodes to handle the switching,
        matching MATLAB's behavior in MNetwork.link().

        References:
            MATLAB: matlab/src/lang/@MNetwork/link.m
        """
        if self._routing_matrix is None:
            return

        from .nodes import ClassSwitch

        nodes = self.get_nodes()
        classes = self.get_classes()
        nclasses = len(classes)
        nnodes = len(nodes)

        if nclasses == 0 or nnodes == 0:
            return

        # Build class switching matrix for each (src_node, dst_node) pair
        # csnodematrix[i][j][r][s] = probability of class r -> class s when going from node i to j
        csnodematrix = {}
        for i in range(nnodes):
            for j in range(nnodes):
                csnodematrix[(i, j)] = np.zeros((nclasses, nclasses))

        # Populate from routing matrix
        for (class_src, class_dst), routes in self._routing_matrix._routes.items():
            src_class_idx = class_src._index if hasattr(class_src, '_index') else classes.index(class_src)
            dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else classes.index(class_dst)

            for (node_src, node_dst), prob in routes.items():
                src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else nodes.index(node_src)
                dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else nodes.index(node_dst)

                if prob > 0:
                    csnodematrix[(src_node_idx, dst_node_idx)][src_class_idx, dst_class_idx] = prob

        # Normalize each row of the switching matrix and check if non-diagonal
        # (jobs that go from node i to node j, conditional on going to j)
        FINE_TOL = 1e-10
        cs_nodes_to_create = {}  # (src_idx, dst_idx) -> normalized switching matrix

        for (i, j), csmat in csnodematrix.items():
            # Normalize rows
            for r in range(nclasses):
                row_sum = np.sum(csmat[r, :])
                if row_sum > FINE_TOL:
                    csmat[r, :] = csmat[r, :] / row_sum
                else:
                    csmat[r, r] = 1.0  # Default: stay in same class

            # Check if non-diagonal (has actual class switching)
            is_diagonal = True
            for r in range(nclasses):
                for s in range(nclasses):
                    if r != s and csmat[r, s] > FINE_TOL:
                        is_diagonal = False
                        break
                if not is_diagonal:
                    break

            if not is_diagonal:
                cs_nodes_to_create[(i, j)] = csmat.copy()

        if not cs_nodes_to_create:
            return  # No class switching needed

        # Create ClassSwitch nodes
        cs_node_map = {}  # (src_idx, dst_idx) -> ClassSwitch node
        for (src_idx, dst_idx), csmat in cs_nodes_to_create.items():
            src_name = nodes[src_idx].name
            dst_name = nodes[dst_idx].name
            cs_name = f'CS_{src_name}_to_{dst_name}'
            cs_node = ClassSwitch(self, cs_name, csmat)
            cs_node._auto_added = True  # Mark as auto-generated
            cs_node_map[(src_idx, dst_idx)] = cs_node

        # Update routing matrix to go through CS nodes
        # Original: P[r,s](i,j) = prob
        # After: P[r,r](i,CS) = sum_s P[r,s](i,j), P[s,s](CS,j) = 1
        new_routes = {}

        for (class_src, class_dst), routes in self._routing_matrix._routes.items():
            src_class_idx = class_src._index if hasattr(class_src, '_index') else classes.index(class_src)
            dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else classes.index(class_dst)

            for (node_src, node_dst), prob in routes.items():
                src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else nodes.index(node_src)
                dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else nodes.index(node_dst)

                if prob > 0 and (src_node_idx, dst_node_idx) in cs_node_map:
                    # Route through CS node
                    cs_node = cs_node_map[(src_node_idx, dst_node_idx)]

                    # P[r,r](src, CS) += P[r,s](src, dst)
                    key = (class_src, class_src)
                    if key not in new_routes:
                        new_routes[key] = {}
                    route_key = (node_src, cs_node)
                    new_routes[key][route_key] = new_routes[key].get(route_key, 0) + prob

                    # P[s,s](CS, dst) = 1 for all classes that can exit CS
                    key = (class_dst, class_dst)
                    if key not in new_routes:
                        new_routes[key] = {}
                    route_key = (cs_node, node_dst)
                    new_routes[key][route_key] = 1.0
                else:
                    # Keep original route (no class switching on this edge)
                    key = (class_src, class_dst)
                    if key not in new_routes:
                        new_routes[key] = {}
                    route_key = (node_src, node_dst)
                    new_routes[key][route_key] = new_routes[key].get(route_key, 0) + prob

        # Update routing matrix
        self._routing_matrix._routes = new_routes
        self._routing_matrix._matrix = None  # Clear cached matrix

    def get_routing_matrix(self) -> Optional[RoutingMatrix]:
        """
        Get the current routing matrix.

        Returns:
            RoutingMatrix or None if not set
        """
        return self._routing_matrix

    def get_linked_routing_matrix(self):
        """
        Get the linked routing matrix (rtorig from NetworkStruct).

        This returns the original routing matrix in the format used
        by the compiled NetworkStruct. The matrix is structured as
        a list of lists: P[r][s] is the routing matrix from class r to class s.

        Returns:
            List of routing matrices, or None if not available
        """
        if not self._has_struct:
            self.refresh_struct(False)

        if self._sn is not None and hasattr(self._sn, 'rtorig') and self._sn.rtorig is not None:
            return self._sn.rtorig

        # Build rtorig from routing matrix if not available
        if self._routing_matrix is None:
            return None

        nclasses = len(self._classes)
        nnodes = len(self._nodes)

        # Create rtorig structure: list of lists of (nnodes x nnodes) matrices
        rtorig = [[np.zeros((nnodes, nnodes)) for _ in range(nclasses)] for _ in range(nclasses)]

        for (from_class, to_class), routes in self._routing_matrix._routes.items():
            from_idx = from_class._index if hasattr(from_class, '_index') else self._classes.index(from_class)
            to_idx = to_class._index if hasattr(to_class, '_index') else self._classes.index(to_class)

            for (from_node, to_node), prob in routes.items():
                from_node_idx = from_node._node_index if hasattr(from_node, '_node_index') else self._nodes.index(from_node)
                to_node_idx = to_node._node_index if hasattr(to_node, '_node_index') else self._nodes.index(to_node)
                rtorig[from_idx][to_idx][from_node_idx, to_node_idx] = prob

        return rtorig

    def reset_network(self, hard: bool = True) -> None:
        """
        Reset the network configuration.

        This clears the routing matrix and struct, allowing the network
        to be reconfigured. Used by model transformations.

        Args:
            hard: If True, also clears the routing matrix
        """
        self._has_struct = False
        self._sn = None
        self._connections = None

        if hard:
            self._routing_matrix = None

    def _refresh_routing_matrix(self) -> None:
        """
        Refresh and validate the routing matrix.

        This method validates that routing probabilities are valid
        (non-negative, stochastic).
        """
        if self._routing_matrix is None:
            return

        # TODO: Add validation for stochasticity, ergodicity
        # This would call dtmc_stochcomp from MATLAB

    def get_connection_matrix(self) -> np.ndarray:
        """
        Get the network connection (adjacency) matrix.

        Returns:
            (N x N) binary matrix where 1 indicates a connection
        """
        if self._connections is None:
            self._compute_connections()
        return self._connections

    def _compute_connections(self) -> None:
        """Compute adjacency matrix from routing matrix."""
        nnodes = len(self._nodes)
        connections = np.zeros((nnodes, nnodes), dtype=int)

        if self._routing_matrix is not None:
            routes = self._routing_matrix.get_routes()
            for (from_class, to_class), route in routes.items():
                # If any route exists from i to j, create edge
                connections = np.logical_or(connections, (route > 0).astype(int))

        self._connections = connections

    # =====================================================================
    # QUERY METHODS
    # =====================================================================

    def get_nodes(self) -> List[Node]:
        """Get list of all nodes."""
        return self._nodes.copy()

    def nodes(self) -> List[Node]:
        """Get list of all nodes (method alias for compatibility)."""
        return self._nodes.copy()

    def get_stations(self) -> List[Station]:
        """Get list of all station nodes."""
        return self._stations.copy()

    def get_classes(self) -> List[JobClass]:
        """Get list of all job classes."""
        return self._classes.copy()

    def classes(self) -> List[JobClass]:
        """Get list of all job classes (method alias for compatibility)."""
        return self._classes.copy()

    def get_number_of_nodes(self) -> int:
        """Get total number of nodes."""
        return len(self._nodes)

    def get_number_of_stations(self) -> int:
        """Get number of station nodes."""
        return len(self._stations)

    def get_number_of_classes(self) -> int:
        """Get number of job classes."""
        return len(self._classes)

    def get_number_of_jobs(self) -> np.ndarray:
        """
        Get population vector for all classes.

        Returns:
            (K,) array with population for each class
        """
        njobs = np.zeros(len(self._classes))
        for i, jobclass in enumerate(self._classes):
            if jobclass.jobclass_type == JobClassType.CLOSED:
                # TODO: Get population from ClosedClass
                pass
        return njobs

    def get_node_by_name(self, name: str) -> Optional[Node]:
        """
        Get node by name.

        Args:
            name: Node name

        Returns:
            Node or None if not found
        """
        for node in self._nodes:
            if node.name == name:
                return node
        return None

    def get_station_by_name(self, name: str) -> Optional[Station]:
        """
        Get station by name.

        Args:
            name: Station name

        Returns:
            Station or None if not found
        """
        for station in self._stations:
            if station.name == name:
                return station
        return None

    def get_class_by_name(self, name: str) -> Optional[JobClass]:
        """
        Get job class by name.

        Args:
            name: Class name

        Returns:
            JobClass or None if not found
        """
        for jobclass in self._classes:
            if jobclass.name == name:
                return jobclass
        return None

    def get_node_index(self, node: Union[Node, str]) -> int:
        """
        Get 1-based index of a node.

        Args:
            node: Node instance or node name

        Returns:
            1-based node index
        """
        if isinstance(node, str):
            node = self.get_node_by_name(node)
            if node is None:
                raise ValueError(f"[{self.name}] Node '{node}' not found")

        try:
            idx = self._nodes.index(node)
            return idx + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(f"[{self.name}] Node '{node.name}' not in network")

    def get_station_index(self, station: Union[Station, str]) -> int:
        """
        Get 1-based index of a station.

        Args:
            station: Station instance or station name

        Returns:
            1-based station index
        """
        if isinstance(station, str):
            station = self.get_station_by_name(station)
            if station is None:
                raise ValueError(f"[{self.name}] Station '{station}' not found")

        try:
            idx = self._stations.index(station)
            return idx + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(f"[{self.name}] Station '{station.name}' not in network")

    def get_class_index(self, jobclass: Union[JobClass, str]) -> int:
        """
        Get 1-based index of a job class.

        Args:
            jobclass: JobClass instance or class name

        Returns:
            1-based class index
        """
        if isinstance(jobclass, str):
            jobclass = self.get_class_by_name(jobclass)
            if jobclass is None:
                raise ValueError(f"[{self.name}] Class '{jobclass}' not found")

        try:
            idx = self._classes.index(jobclass)
            return idx + 1  # Convert to 1-based
        except ValueError:
            raise ValueError(f"[{self.name}] Class '{jobclass.name}' not in network")

    def get_source(self) -> Optional[Node]:
        """Get Source node (if present)."""
        if self._source_idx >= 0:
            return self._nodes[self._source_idx]

        for node in self._nodes:
            if node.node_type == NodeType.SOURCE:
                self._source_idx = self._nodes.index(node)
                return node

        return None

    def get_sink(self) -> Optional[Node]:
        """Get Sink node (if present)."""
        if self._sink_idx >= 0:
            return self._nodes[self._sink_idx]

        for node in self._nodes:
            if node.node_type == NodeType.SINK:
                self._sink_idx = self._nodes.index(node)
                return node

        return None

    def has_open_classes(self) -> bool:
        """Check if network has any open classes."""
        return any(cls.jobclass_type == JobClassType.OPEN for cls in self._classes)

    def has_closed_classes(self) -> bool:
        """Check if network has any closed classes."""
        return any(cls.jobclass_type == JobClassType.CLOSED for cls in self._classes)

    def has_fork(self) -> bool:
        """Check if network has any fork nodes."""
        from .nodes import Fork
        return any(isinstance(node, Fork) for node in self._nodes)

    def has_join(self) -> bool:
        """Check if network has any join nodes."""
        from .nodes import Join
        return any(isinstance(node, Join) for node in self._nodes)

    # =====================================================================
    # STRUCT COMPILATION METHODS
    # =====================================================================

    def refresh_struct(self, hard: bool = True) -> None:
        """
        Compile model to NetworkStruct for solver consumption.

        Args:
            hard: If True, recompute everything; if False, update only changed parts
        """
        if self._sn is None:
            hard = True

        if hard:
            self._build_struct_from_scratch()
        else:
            # TODO: Implement incremental updates
            self._build_struct_from_scratch()

        self._has_struct = True

    def get_struct(self) -> NetworkStruct:
        """
        Get compiled NetworkStruct.

        Returns:
            NetworkStruct for solver use

        Raises:
            ValueError: If model hasn't been compiled
        """
        if not self._has_struct:
            self.refresh_struct()

        return self._sn

    def reset_struct(self) -> None:
        """Reset struct compilation flag."""
        self._reset_struct()

    def reset(self) -> None:
        """
        Reset the network to allow re-configuration.

        This resets the routing matrix and struct, allowing nodes to be
        reconfigured (e.g., changing routing strategies) before re-solving.
        """
        self._has_struct = False
        self._sn = None
        self._connections = None

    def struct(self) -> 'NetworkStruct':
        """
        Get compiled NetworkStruct (alias for get_struct).

        Returns:
            NetworkStruct for solver use
        """
        return self.get_struct()

    def relink(self, routing_matrix: RoutingMatrix) -> None:
        """
        Re-link the network with a new routing matrix.

        This is similar to link() but resets the struct first to allow
        iterative optimization of routing parameters.

        Args:
            routing_matrix: New RoutingMatrix specifying routing probabilities
        """
        self._has_struct = False
        self._sn = None
        self.link(routing_matrix)

    def _reset_struct(self) -> None:
        """Internal method to invalidate struct."""
        self._has_struct = False
        self._sn = None

    def _build_struct_from_scratch(self) -> None:
        """
        Build NetworkStruct from scratch.

        This is the main compilation method that converts the network
        into a form suitable for solvers.
        """
        # Validate network
        self._validate()

        # Initialize struct
        self._sn = NetworkStruct()

        # Set basic dimensions
        self._sn.nstations = len(self._stations)
        self._sn.nclasses = len(self._classes)
        self._sn.nnodes = len(self._nodes)

        # Build node mappings
        self._refresh_node_mappings()

        # Extract service/arrival rates
        self._refresh_rates()

        # Extract load-dependent scaling
        self._refresh_load_dependence()

        # Extract scheduling strategies
        self._refresh_scheduling()

        # Extract capacity
        self._refresh_capacity()

        # Build routing matrix
        self._refresh_routing()

        # Compute chains and visits
        self._refresh_chains()

        # Extract node parameters (for transitions, caches, etc.)
        self._refresh_nodeparam()

        # Build fork-join relationship matrix
        self._refresh_fork_joins()

        # Extract initial state (for SPNs)
        self._refresh_state()

    def _validate(self) -> None:
        """
        Validate network structure.

        Raises:
            ValueError: If network is invalid
        """
        if len(self._nodes) == 0:
            raise ValueError(f"[{self.name}] Network has no nodes")

        if len(self._classes) == 0:
            raise ValueError(f"[{self.name}] Network has no job classes")

        # Check for required nodes (Source/Sink for open classes)
        if self.has_open_classes():
            if self.get_source() is None:
                raise ValueError(f"[{self.name}] Open classes require Source node")
            if self.get_sink() is None:
                raise ValueError(f"[{self.name}] Open classes require Sink node")

    def _refresh_node_mappings(self) -> None:
        """Build mappings between nodes and stations."""
        nnodes = len(self._nodes)
        nstations = len(self._stations)
        nclasses = len(self._classes)

        self._sn.nodenames = [node.name for node in self._nodes]
        self._sn.classnames = [cls.name for cls in self._classes]

        # Node types
        self._sn.nodetype = [node.node_type for node in self._nodes]

        # Station classification
        self._sn.isstation = np.array([node.is_station() for node in self._nodes], dtype=bool)

        # Stateful classification (nodes that maintain state)
        # Stateful: SOURCE, DELAY, QUEUE, CACHE, JOIN, ROUTER, PLACE, TRANSITION
        from .base import NodeType
        stateful_types = {NodeType.SOURCE, NodeType.DELAY, NodeType.QUEUE,
                         NodeType.CACHE, NodeType.JOIN, NodeType.ROUTER,
                         NodeType.PLACE, NodeType.TRANSITION}
        self._sn.isstateful = np.array([node.node_type in stateful_types for node in self._nodes], dtype=bool)
        self._sn.nstateful = int(np.sum(self._sn.isstateful))

        # Node to station mapping
        node_to_station = np.full(nnodes, -1, dtype=int)
        for i, station in enumerate(self._stations):
            node_idx = self._nodes.index(station)
            node_to_station[node_idx] = i
        self._sn.nodeToStation = node_to_station

        # Station to node mapping
        station_to_node = np.zeros(nstations, dtype=int)
        for i, station in enumerate(self._stations):
            station_to_node[i] = self._nodes.index(station)
        self._sn.stationToNode = station_to_node

        # Node to stateful mapping
        node_to_stateful = np.full(nnodes, -1, dtype=int)
        stateful_idx = 0
        for i, node in enumerate(self._nodes):
            if self._sn.isstateful[i]:
                node_to_stateful[i] = stateful_idx
                stateful_idx += 1
        self._sn.nodeToStateful = node_to_stateful

        # Stateful to node mapping
        stateful_to_node = np.zeros(self._sn.nstateful, dtype=int)
        stateful_to_station = np.full(self._sn.nstateful, -1, dtype=int)
        stateful_idx = 0
        for i, node in enumerate(self._nodes):
            if self._sn.isstateful[i]:
                stateful_to_node[stateful_idx] = i
                if node.is_station():
                    stateful_to_station[stateful_idx] = node_to_station[i]
                stateful_idx += 1
        self._sn.statefulToNode = stateful_to_node
        self._sn.statefulToStation = stateful_to_station

        # Station to stateful mapping
        station_to_stateful = np.zeros(nstations, dtype=int)
        for i, station in enumerate(self._stations):
            node_idx = self._nodes.index(station)
            station_to_stateful[i] = node_to_stateful[node_idx]
        self._sn.stationToStateful = station_to_stateful

        # Number of servers per station
        nservers = np.ones(nstations)
        for i, station in enumerate(self._stations):
            if hasattr(station, '_number_of_servers'):
                nservers[i] = station._number_of_servers
            elif hasattr(station, 'get_number_of_servers'):
                nservers[i] = station.get_number_of_servers()
        self._sn.nservers = nservers

        # Population per class
        njobs = np.zeros(nclasses)
        refstat = np.zeros(nclasses, dtype=int)
        classprio = np.zeros(nclasses)

        for j, jobclass in enumerate(self._classes):
            # Get priority
            if hasattr(jobclass, '_priority'):
                classprio[j] = jobclass._priority

            # Get reference station
            ref = jobclass.get_reference_station() if hasattr(jobclass, 'get_reference_station') else None
            if ref is None:
                ref = getattr(jobclass, '_refstat', None)
            if ref is None:
                ref = getattr(jobclass, '_reference_station', None)

            if ref is not None and ref in self._stations:
                refstat[j] = self._stations.index(ref)
            else:
                # For open classes, reference is typically source (station 0)
                refstat[j] = 0

            # Get population (closed classes)
            if hasattr(jobclass, '_njobs'):
                njobs[j] = jobclass._njobs
            elif hasattr(jobclass, 'getNumberOfJobs'):
                njobs[j] = jobclass.getNumberOfJobs()
            else:
                # Open class - infinite population
                njobs[j] = np.inf

        self._sn.njobs = njobs
        self._sn.refstat = refstat
        self._sn.classprio = classprio
        self._sn.nclosedjobs = int(np.sum(njobs[np.isfinite(njobs)]))

    def _refresh_rates(self) -> None:
        """Extract service and arrival rates from nodes."""
        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Service rates: (M x K) array - rate = 1/mean_service_time
        self._sn.rates = np.zeros((nstations, nclasses))

        # Squared coefficient of variation
        self._sn.scv = np.ones((nstations, nclasses))

        # Process parameters: nested list [station][class] -> [D0, D1] or similar
        self._sn.proc = [[None for _ in range(nclasses)] for _ in range(nstations)]

        # Process type IDs: (M x K) array of ProcessType values
        self._sn.procid = np.empty((nstations, nclasses), dtype=object)
        self._sn.procid.fill(ProcessType.EXP)  # Default to exponential

        for i, station in enumerate(self._stations):
            for j, jobclass in enumerate(self._classes):
                dist = None

                # Try to get service distribution
                if hasattr(station, 'get_service'):
                    dist = station.get_service(jobclass)
                elif hasattr(station, '_service_process'):
                    dist = station._service_process.get(jobclass)

                # For Source nodes, get arrival distribution
                if dist is None and hasattr(station, 'get_arrival'):
                    dist = station.get_arrival(jobclass)
                if dist is None and hasattr(station, '_arrival_process'):
                    dist = getattr(station, '_arrival_process', {}).get(jobclass)

                if dist is not None:
                    # Extract mean service time
                    mean = None
                    if hasattr(dist, 'getMean'):
                        mean = dist.getMean()
                    elif hasattr(dist, 'mean'):
                        mean = dist.mean
                    elif hasattr(dist, '_mean'):
                        mean = dist._mean

                    if mean is not None and mean > 0:
                        self._sn.rates[i, j] = 1.0 / mean

                    # Extract squared coefficient of variation
                    scv = None
                    if hasattr(dist, 'getSCV'):
                        scv = dist.getSCV()
                    elif hasattr(dist, 'scv'):
                        scv = dist.scv
                    elif hasattr(dist, '_scv'):
                        scv = dist._scv

                    if scv is not None:
                        self._sn.scv[i, j] = scv

                    # Extract process type and parameters
                    self._extract_process_params(i, j, dist)
                else:
                    # No distribution defined
                    # Join nodes should have rate=0 (they don't do service work)
                    from .nodes import Join
                    if isinstance(station, Join):
                        self._sn.rates[i, j] = 0.0
                    else:
                        # Default rate of 1.0 for other stations
                        self._sn.rates[i, j] = 1.0

    def _extract_process_params(self, station_idx: int, class_idx: int, dist) -> None:
        """
        Extract process type and parameters from a distribution.

        Populates self._sn.proc and self._sn.procid for the given station/class.

        Args:
            station_idx: Station index
            class_idx: Class index
            dist: Distribution object
        """
        # Import distribution classes for type checking
        from ..distributions.markovian import MAP, MMPP2, PH, APH, Coxian
        from ..distributions.continuous import (
            Exp, Erlang, HyperExp, Gamma, Uniform, Det, Pareto, Weibull, Lognormal
        )

        # Check for MMPP2 first (more specific than MAP)
        if isinstance(dist, MMPP2):
            self._sn.procid[station_idx, class_idx] = ProcessType.MMPP2
            D0 = dist.getD0()
            D1 = dist.getD1()
            self._sn.proc[station_idx][class_idx] = [D0, D1]

        # Check for MAP
        elif isinstance(dist, MAP):
            self._sn.procid[station_idx, class_idx] = ProcessType.MAP
            D0 = dist.getD0()
            D1 = dist.getD1()
            self._sn.proc[station_idx][class_idx] = [D0, D1]

        # Check for Phase-Type distributions
        elif isinstance(dist, (PH, APH)):
            if isinstance(dist, APH):
                self._sn.procid[station_idx, class_idx] = ProcessType.APH
            else:
                self._sn.procid[station_idx, class_idx] = ProcessType.PH
            # Store [alpha, T] where alpha is initial prob vector, T is generator
            alpha = dist.getInitProb() if hasattr(dist, 'getInitProb') else None
            T = dist.getD0() if hasattr(dist, 'getD0') else None
            if alpha is not None and T is not None:
                self._sn.proc[station_idx][class_idx] = [alpha, T]

        elif isinstance(dist, Coxian):
            self._sn.procid[station_idx, class_idx] = ProcessType.COXIAN
            alpha = dist.getInitProb() if hasattr(dist, 'getInitProb') else None
            T = dist.getD0() if hasattr(dist, 'getD0') else None
            if alpha is not None and T is not None:
                self._sn.proc[station_idx][class_idx] = [alpha, T]

        # Check for HyperExp
        elif isinstance(dist, HyperExp):
            self._sn.procid[station_idx, class_idx] = ProcessType.HYPEREXP
            # Store parameters for HyperExp (uses _probs and _rates arrays)
            if hasattr(dist, '_probs') and hasattr(dist, '_rates'):
                self._sn.proc[station_idx][class_idx] = {
                    'probs': dist._probs, 'rates': dist._rates
                }

        # Check for Erlang
        elif isinstance(dist, Erlang):
            self._sn.procid[station_idx, class_idx] = ProcessType.ERLANG
            # Erlang uses _phases and _phase_rate attributes
            if hasattr(dist, '_phases') and hasattr(dist, '_phase_rate'):
                self._sn.proc[station_idx][class_idx] = {
                    'k': dist._phases, 'mu': dist._phase_rate
                }

        # Check for Exponential
        elif isinstance(dist, Exp):
            self._sn.procid[station_idx, class_idx] = ProcessType.EXP
            # For Exp, proc can remain None or store rate
            if hasattr(dist, '_rate'):
                self._sn.proc[station_idx][class_idx] = {'rate': dist._rate}

        # Check for other distributions
        elif isinstance(dist, Gamma):
            self._sn.procid[station_idx, class_idx] = ProcessType.GAMMA
        elif isinstance(dist, Uniform):
            self._sn.procid[station_idx, class_idx] = ProcessType.UNIFORM
        elif isinstance(dist, Det):
            self._sn.procid[station_idx, class_idx] = ProcessType.DET
        elif isinstance(dist, Pareto):
            self._sn.procid[station_idx, class_idx] = ProcessType.PARETO
        elif isinstance(dist, Weibull):
            self._sn.procid[station_idx, class_idx] = ProcessType.WEIBULL
        elif isinstance(dist, Lognormal):
            self._sn.procid[station_idx, class_idx] = ProcessType.LOGNORMAL

        # Fallback: try to detect by class name
        else:
            class_name = type(dist).__name__
            proc_type = ProcessType.fromString(class_name)
            if proc_type is not None:
                self._sn.procid[station_idx, class_idx] = proc_type
            # Keep default EXP if unknown

    def _normalize_sched_strategy(self, sched) -> SchedStrategy:
        """
        Normalize scheduling strategy from any SchedStrategy enum to native format.

        Different modules define SchedStrategy with different integer values.
        This normalizes by matching on the strategy name.

        Args:
            sched: SchedStrategy enum from any source

        Returns:
            SchedStrategy from lang.base with correct native value
        """
        # Get the name of the scheduling strategy
        if hasattr(sched, 'name'):
            name = sched.name
        else:
            name = str(sched).split('.')[-1]

        # Map to native SchedStrategy by name
        try:
            return SchedStrategy[name]
        except KeyError:
            # Fallback to FCFS if unknown
            return SchedStrategy.FCFS

    def _refresh_load_dependence(self) -> None:
        """Extract load-dependent scaling from stations."""
        nstations = len(self._stations)

        # Get total population for closed classes
        njobs = self._sn.njobs
        total_pop = int(np.sum(njobs[np.isfinite(njobs)]))

        if total_pop == 0:
            # No closed jobs, no load dependence needed
            self._sn.lldscaling = None
            return

        # Check if any station has load dependence
        has_ld = False
        for station in self._stations:
            ld = station.get_load_dependence() if hasattr(station, 'get_load_dependence') else None
            if ld is not None and len(ld) > 0:
                has_ld = True
                break

        if not has_ld:
            # No load dependence configured
            self._sn.lldscaling = None
            return

        # Build lldscaling matrix: (M x Nmax)
        # lldscaling[station, n-1] = scaling factor when n jobs are at that station
        lldscaling = np.ones((nstations, total_pop))

        for i, station in enumerate(self._stations):
            ld = station.get_load_dependence() if hasattr(station, 'get_load_dependence') else None
            if ld is not None and len(ld) > 0:
                ld_array = np.asarray(ld)
                # Copy scaling factors, repeating last value if needed
                for j in range(total_pop):
                    if j < len(ld_array):
                        lldscaling[i, j] = ld_array[j]
                    else:
                        lldscaling[i, j] = ld_array[-1]

        self._sn.lldscaling = lldscaling

    def _refresh_scheduling(self) -> None:
        """Extract scheduling strategies from stations."""
        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Scheduling strategies
        self._sn.sched = {}
        for i, station in enumerate(self._stations):
            sched = SchedStrategy.FCFS  # Default

            # Try to get scheduling strategy
            if hasattr(station, 'get_sched_strategy'):
                sched = station.get_sched_strategy()
            elif hasattr(station, '_sched_strategy'):
                sched = station._sched_strategy
            elif hasattr(station, 'node_type'):
                # Delay nodes use INF scheduling
                from .base import NodeType
                if station.node_type == NodeType.DELAY:
                    sched = SchedStrategy.INF
                elif station.node_type == NodeType.SOURCE:
                    sched = SchedStrategy.EXT

            # Normalize to native format (different modules use different enum values)
            sched = self._normalize_sched_strategy(sched)
            self._sn.sched[i] = sched

        # Scheduling parameters (weights, priorities for DPS/GPS/PS)
        self._sn.schedparam = np.ones((nstations, nclasses))

        for i, station in enumerate(self._stations):
            for j, jobclass in enumerate(self._classes):
                param = 1.0  # Default weight

                if hasattr(station, 'get_strategy_param'):
                    param = station.get_strategy_param(jobclass)
                elif hasattr(station, '_sched_param'):
                    param = station._sched_param.get(jobclass, 1.0)

                self._sn.schedparam[i, j] = param

    def _refresh_capacity(self) -> None:
        """Extract capacity limits from stations."""
        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Total capacity per station
        self._sn.cap = np.array([station.capacity for station in self._stations])

        # Per-class capacity
        classcap = np.full((nstations, nclasses), np.inf)
        for i, station in enumerate(self._stations):
            for j, jobclass in enumerate(self._classes):
                classcap[i, j] = station.get_class_capacity(jobclass)

        self._sn.classcap = classcap

    def _refresh_routing(self) -> None:
        """Build routing matrix in NetworkStruct."""
        nstations = len(self._stations)
        nclasses = len(self._classes)
        nnodes = len(self._nodes)

        if nstations == 0 or nclasses == 0:
            return

        # Get station-indexed routing matrix from RoutingMatrix
        if self._routing_matrix is not None:
            # toMatrix() returns (M*K) x (M*K) station-indexed matrix
            self._sn.rt = self._routing_matrix.toMatrix()
        else:
            # No routing defined - create empty matrix
            self._sn.rt = np.zeros((nstations * nclasses, nstations * nclasses))

        # Build node-indexed routing matrix (rtnodes)
        # This includes all nodes, not just stations
        self._sn.rtnodes = np.zeros((nnodes * nclasses, nnodes * nclasses))

        if self._routing_matrix is not None:
            # Copy from sparse routes to node-indexed matrix
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                src_class_idx = class_src._index if hasattr(class_src, '_index') else self._classes.index(class_src)
                dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else self._classes.index(class_dst)

                for (node_src, node_dst), prob in routes.items():
                    src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                    dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)

                    i = src_node_idx * nclasses + src_class_idx
                    j = dst_node_idx * nclasses + dst_class_idx
                    self._sn.rtnodes[i, j] = prob

        # Apply ClassSwitch matrices to rtnodes
        # This matches MATLAB's getRoutingMatrix behavior for StatelessClassSwitcher
        from .nodes import ClassSwitch
        for ind, node in enumerate(self._nodes):
            if isinstance(node, ClassSwitch) and hasattr(node, '_switch_matrix') and node._switch_matrix is not None:
                Pcs = np.asarray(node._switch_matrix)
                if Pcs.shape[0] == nclasses and Pcs.shape[1] == nclasses:
                    # Extract routing rows for this node (make a copy to avoid modification)
                    Pi = self._sn.rtnodes[ind * nclasses:(ind + 1) * nclasses, :].copy()

                    # Zero out original routing from this node
                    self._sn.rtnodes[ind * nclasses:(ind + 1) * nclasses, :] = 0

                    # Apply class switching: for each destination node jnd
                    for jnd in range(nnodes):
                        # Pij[r,s] = Pi[r, (jnd-1)*K+s] - routing from class r at ind to class s at jnd
                        Pij = Pi[:, jnd * nclasses:(jnd + 1) * nclasses]
                        # diag(Pij) gives routing probabilities for same-class transitions
                        diag_Pij = np.diag(Pij)
                        # New routing: Pcs[r,s] * Pij[s,s] for all (r,s)
                        # This means: class r switches to class s (Pcs[r,s]) and then routes as class s (Pij[s,s])
                        new_routing = Pcs * diag_Pij[np.newaxis, :]  # broadcast diag across rows
                        self._sn.rtnodes[ind * nclasses:(ind + 1) * nclasses,
                                        jnd * nclasses:(jnd + 1) * nclasses] = new_routing

        # Compute station-indexed routing matrix (rt) from node-indexed (rtnodes)
        # using stochastic complement to absorb non-station nodes
        # This matches MATLAB: rt = dtmc_stochcomp(rtnodes, statefulNodesClasses)
        from ..api.mc.dtmc import dtmc_stochcomp

        # Build list of state indices corresponding to stations (stateful nodes)
        stateful_nodes_classes = []
        for node in self._nodes:
            if hasattr(node, '_station_index') and node._station_index is not None and node._station_index >= 0:
                node_idx = node._node_index if hasattr(node, '_node_index') else self._nodes.index(node)
                for k in range(nclasses):
                    stateful_nodes_classes.append(node_idx * nclasses + k)

        if len(stateful_nodes_classes) > 0:
            stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)
            # Compute stochastic complement
            rt_from_rtnodes = dtmc_stochcomp(self._sn.rtnodes, stateful_nodes_classes)
            # Reorder from node-class indexing to station-class indexing
            # The stochcomp result is indexed by the keep_states order, which is
            # station nodes in node order (which should match station order)
            self._sn.rt = rt_from_rtnodes

        # Build connection matrix
        self._sn.connmatrix = np.zeros((nnodes, nnodes))
        if self._routing_matrix is not None:
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                for (node_src, node_dst), prob in routes.items():
                    if prob > 0:
                        src_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                        dst_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)
                        self._sn.connmatrix[src_idx, dst_idx] = 1

        # Build routing strategy matrix (N x K) - stores routing strategy per node/class
        from .base import RoutingStrategy
        self._sn.routing = np.full((nnodes, nclasses), int(RoutingStrategy.RAND), dtype=int)
        for i, node in enumerate(self._nodes):
            if hasattr(node, '_routing_strategies') and node._routing_strategies:
                for jobclass, strategy in node._routing_strategies.items():
                    class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                    # Handle both IntEnum and plain int values
                    if hasattr(strategy, 'value'):
                        self._sn.routing[i, class_idx] = strategy.value
                    else:
                        self._sn.routing[i, class_idx] = int(strategy)

    def _refresh_chains(self) -> None:
        """
        Compute chains and visit ratios.

        A chain is a group of job classes that can transition into each other
        via class switching. Classes in the same chain share routing structure.
        """
        nstations = len(self._stations)
        nclasses = len(self._classes)

        if nstations == 0 or nclasses == 0:
            self._sn.nchains = 0
            self._sn.chains = np.array([])
            self._sn.visits = {}
            self._sn.inchain = {}
            return

        # Build class transition graph (csmask) to detect chains
        # Two classes are in the same chain if there's routing from one to the other.
        # We use the rt matrix which already incorporates class switching through
        # ClassSwitch nodes via stochastic complement.
        # This matches MATLAB's getRoutingMatrix.m which builds csmask from rt.
        class_connected = np.zeros((nclasses, nclasses), dtype=bool)

        if self._sn.rt is not None and self._sn.rt.size > 0:
            rt = self._sn.rt
            M = nstations
            K = nclasses
            # csmask[r,s] = True if any station routes class r to class s
            for r in range(K):
                for s in range(K):
                    for ist in range(M):
                        for jst in range(M):
                            if rt[ist * K + r, jst * K + s] > 0:
                                class_connected[r, s] = True
                                break
                        if class_connected[r, s]:
                            break

        # Find connected components using union-find
        parent = list(range(nclasses))

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        # Union classes that can transition to each other
        # This forms chains based on class switching - classes that can route
        # to each other (directly or indirectly) are in the same chain.
        # NOTE: We do NOT union classes based on reference station, as MATLAB
        # uses weaklyconncomp on csmask to determine chains purely from class
        # switching structure (see getRoutingMatrix.m lines 279-291).
        for i in range(nclasses):
            for j in range(nclasses):
                if class_connected[i, j] or class_connected[j, i]:
                    union(i, j)

        # Group classes by chain
        chain_classes = {}
        for i in range(nclasses):
            root = find(i)
            if root not in chain_classes:
                chain_classes[root] = []
            chain_classes[root].append(i)

        # Assign chain IDs
        chain_roots = sorted(chain_classes.keys())
        nchains = len(chain_roots)
        chains = np.zeros(nclasses, dtype=int)
        inchain = {}

        for chain_id, root in enumerate(chain_roots):
            class_indices = chain_classes[root]
            for class_idx in class_indices:
                chains[class_idx] = chain_id
            inchain[chain_id] = np.array(class_indices, dtype=int)

        self._sn.nchains = nchains
        self._sn.chains = chains
        self._sn.inchain = inchain

        # Compute visit ratios for each chain
        self._sn.visits = {}

        if self._sn.rt is not None and self._sn.rt.size > 0:
            from ..api.mc import dtmc_solve

            for chain_id in range(nchains):
                classes_in_chain = inchain[chain_id]
                n_chain_classes = len(classes_in_chain)

                # Check if this is an open chain (has infinite population)
                is_open_chain = False
                njobs = self._sn.njobs.flatten() if self._sn.njobs is not None else np.zeros(nclasses)
                for k in classes_in_chain:
                    if k < len(njobs) and np.isinf(njobs[k]):
                        is_open_chain = True
                        break

                # Build submatrix for this chain
                # Extract rows/cols for (station, class) pairs in this chain
                indices = []
                for ist in range(nstations):
                    for k in classes_in_chain:
                        indices.append(ist * nclasses + k)

                if len(indices) > 0:
                    P_chain = self._sn.rt[np.ix_(indices, indices)]

                    if is_open_chain:
                        # For open chains, compute visits using traffic equations
                        # v = v * P + lambda, where lambda has 1 at Source station
                        # This gives v = lambda * (I - P)^{-1}
                        # For simplicity, we propagate visits from Source through the network
                        n = len(indices)
                        visits_vec = np.zeros(n)

                        # Find Source station (EXT scheduling)
                        source_idx = -1
                        if self._sn.sched is not None:
                            for ist, sched_val in self._sn.sched.items():
                                from ..api.sn import SchedStrategy
                                if sched_val == SchedStrategy.EXT or (isinstance(sched_val, int) and sched_val == 16):
                                    # Map station to index in P_chain
                                    for k in classes_in_chain:
                                        source_idx = ist * n_chain_classes + classes_in_chain.tolist().index(k) if k in classes_in_chain else -1
                                        if source_idx >= 0 and source_idx < n:
                                            break
                                    break

                        if source_idx >= 0 and source_idx < n:
                            # Initialize with arrival at source
                            visits_vec[source_idx] = 1.0

                            # Solve traffic equations: v = e_source * (I - P)^{-1}
                            # where e_source is unit vector at source
                            try:
                                I_minus_P = np.eye(n) - P_chain
                                # Use pseudo-inverse for numerical stability
                                visits_vec = np.linalg.solve(I_minus_P.T, visits_vec)
                                visits_vec = np.abs(visits_vec)  # Ensure non-negative
                            except np.linalg.LinAlgError:
                                # Fallback: propagate through routing iteratively
                                visits_vec = np.zeros(n)
                                visits_vec[source_idx] = 1.0
                                for _ in range(n * 2):
                                    visits_vec = visits_vec + visits_vec @ P_chain
                        else:
                            # No source found, use uniform visits
                            visits_vec = np.ones(n)

                        # Reshape to visits matrix
                        visits_chain = np.zeros((nstations, nclasses))
                        idx = 0
                        for ist in range(nstations):
                            for k_idx, k in enumerate(classes_in_chain):
                                if idx < len(visits_vec):
                                    visits_chain[ist, k] = visits_vec[idx]
                                idx += 1

                        # Normalize by reference station
                        ref_class = classes_in_chain[0]
                        ref_stat = int(self._sn.refstat[ref_class])
                        ref_visit = visits_chain[ref_stat, ref_class]
                        if ref_visit > 0:
                            visits_chain = visits_chain / ref_visit

                        self._sn.visits[chain_id] = visits_chain
                    else:
                        # Closed chain: compute visits
                        n = len(indices)

                        # Check if network is reducible (has absorbing states)
                        has_absorbing = np.any(np.abs(np.diag(P_chain) - 1.0) < 1e-10)

                        if has_absorbing:
                            # Reducible network: use uniform visits (matches MATLAB behavior)
                            # MATLAB assigns equal visits to all stations in reducible networks
                            pi = np.ones(n)
                        else:
                            # Ergodic network: use steady-state distribution
                            row_sums = P_chain.sum(axis=1)
                            visited = row_sums > 1e-10

                            if np.sum(visited) > 0:
                                P_visited = P_chain[np.ix_(np.where(visited)[0], np.where(visited)[0])]
                                row_sums_visited = P_visited.sum(axis=1, keepdims=True)
                                nonzero = (row_sums_visited > 0).flatten()
                                if np.any(nonzero):
                                    P_visited[nonzero] = P_visited[nonzero] / row_sums_visited[nonzero]

                                try:
                                    pi_visited = dtmc_solve(P_visited)
                                except Exception:
                                    pi_visited = np.ones(np.sum(visited)) / np.sum(visited)

                                pi = np.zeros(n)
                                pi[visited] = pi_visited
                            else:
                                pi = np.ones(n) / n

                        # Reshape to (nstations, n_chain_classes)
                        visits_chain = np.zeros((nstations, nclasses))
                        idx = 0
                        for ist in range(nstations):
                            for k_idx, k in enumerate(classes_in_chain):
                                if idx < len(pi):
                                    visits_chain[ist, k] = pi[idx]
                                idx += 1

                        # Normalize by SUM of visits at reference station for ALL classes in chain
                        # This matches MATLAB's sn_refresh_visits.m line 118-121:
                        # normSum = sum(visits{c}(refstat(inchain{c}(1)),inchain{c}))
                        ref_class = classes_in_chain[0]
                        ref_stat = int(self._sn.refstat[ref_class])
                        norm_sum = np.sum(visits_chain[ref_stat, classes_in_chain])
                        if norm_sum > 1e-10:
                            visits_chain = visits_chain / norm_sum

                        self._sn.visits[chain_id] = visits_chain
                else:
                    self._sn.visits[chain_id] = np.ones((nstations, nclasses))
        else:
            # No routing - default visits
            for chain_id in range(nchains):
                self._sn.visits[chain_id] = np.ones((nstations, nclasses))

        # Set reference class per chain
        refclass = np.zeros(nchains, dtype=int)
        for chain_id in range(nchains):
            refclass[chain_id] = inchain[chain_id][0]
        self._sn.refclass = refclass

    def _refresh_nodeparam(self) -> None:
        """
        Extract node parameters for transitions, caches, and other special nodes.

        This populates sn.nodeparam with mode information for transitions in SPNs.
        """
        from .nodes import Transition
        from dataclasses import dataclass, field
        from typing import Optional, List, Dict, Any

        nnodes = len(self._nodes)
        nclasses = len(self._classes)

        # Initialize nodeparam dictionary
        nodeparam = {}

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Transition):
                # Extract transition mode information
                nmodes = node.get_number_of_modes()

                @dataclass
                class TransitionParam:
                    """Container for transition parameters."""
                    nmodes: int = 1
                    modenames: List[str] = field(default_factory=list)
                    timingstrategies: List[Any] = field(default_factory=list)
                    firingprio: List[int] = field(default_factory=list)
                    fireweight: List[float] = field(default_factory=list)
                    nmodeservers: np.ndarray = field(default_factory=lambda: np.array([1.0]))
                    enabling: List[np.ndarray] = field(default_factory=list)
                    inhibiting: List[np.ndarray] = field(default_factory=list)
                    firing: List[np.ndarray] = field(default_factory=list)
                    distributions: List[Any] = field(default_factory=list)

                param = TransitionParam(nmodes=nmodes)

                # Mode names
                param.modenames = node._mode_names.copy() if node._mode_names else [f'Mode{i}' for i in range(nmodes)]

                # Timing strategies
                param.timingstrategies = node._timing_strategies.copy() if node._timing_strategies else ['TIMED'] * nmodes

                # Priorities
                param.firingprio = [int(p) for p in node._firing_priorities] if node._firing_priorities else [0] * nmodes

                # Weights
                param.fireweight = [float(w) for w in node._firing_weights] if node._firing_weights else [1.0] * nmodes

                # Number of servers per mode
                param.nmodeservers = np.array(node._number_of_servers, dtype=float) if node._number_of_servers else np.ones(nmodes)

                # Enabling conditions (list of matrices, one per mode)
                param.enabling = []
                for mode_idx in range(nmodes):
                    if node._enabling_conditions and mode_idx < len(node._enabling_conditions):
                        param.enabling.append(node._enabling_conditions[mode_idx])
                    else:
                        param.enabling.append(np.zeros((nnodes, nclasses)))

                # Inhibiting conditions
                param.inhibiting = []
                for mode_idx in range(nmodes):
                    if node._inhibiting_conditions and mode_idx < len(node._inhibiting_conditions):
                        param.inhibiting.append(node._inhibiting_conditions[mode_idx])
                    else:
                        param.inhibiting.append(np.full((nnodes, nclasses), np.inf))

                # Firing outcomes
                param.firing = []
                for mode_idx in range(nmodes):
                    if node._firing_outcomes and mode_idx < len(node._firing_outcomes):
                        param.firing.append(node._firing_outcomes[mode_idx])
                    else:
                        param.firing.append(np.zeros((nnodes, nclasses)))

                # Distributions
                param.distributions = node._distributions.copy() if node._distributions else [None] * nmodes

                nodeparam[node_idx] = param

        # Store in NetworkStruct
        self._sn.nodeparam = nodeparam if nodeparam else None

    def _refresh_fork_joins(self) -> None:
        """
        Build fork-join relationship matrix.

        The fj matrix is a (nnodes x nnodes) boolean matrix where
        fj[i, j] is True if node j is a Join that synchronizes
        jobs forked by node i (a Fork).
        """
        from .nodes import Fork, Join

        nnodes = len(self._nodes)
        fj = np.zeros((nnodes, nnodes), dtype=bool)

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Join):
                fork = node.get_fork() if hasattr(node, 'get_fork') else node._fork
                if fork is not None and fork in self._nodes:
                    fork_idx = self._nodes.index(fork)
                    fj[fork_idx, node_idx] = True

        self._sn.fj = fj

    def _refresh_state(self) -> None:
        """
        Extract initial state from stateful nodes.

        This populates sn.state with initial token counts for Places in SPNs.
        """
        from .nodes import Place, StatefulNode

        nstateful = self._sn.nstateful if hasattr(self._sn, 'nstateful') else 0
        nclasses = len(self._classes)

        if nstateful == 0:
            self._sn.state = {}
            return

        # Initialize state array
        state = []
        for i in range(nstateful):
            state.append(np.zeros(nclasses))

        # Get state from stateful nodes
        nodeToStateful = self._sn.nodeToStateful if hasattr(self._sn, 'nodeToStateful') else None
        if nodeToStateful is not None:
            nodeToStateful = np.asarray(nodeToStateful).flatten()

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, StatefulNode):
                # Get stateful index
                stateful_idx = None
                if nodeToStateful is not None and node_idx < len(nodeToStateful):
                    stateful_idx = int(nodeToStateful[node_idx])

                if stateful_idx is not None and stateful_idx >= 0 and stateful_idx < len(state):
                    # Get node state
                    node_state = node.get_state() if hasattr(node, 'get_state') else node._state if hasattr(node, '_state') else None
                    if node_state is not None:
                        node_state = np.asarray(node_state).flatten()
                        for k in range(min(nclasses, len(node_state))):
                            state[stateful_idx][k] = node_state[k]

        self._sn.state = state

    # =====================================================================
    # REWARD METHODS
    # =====================================================================

    def set_reward(self, name: str, reward_fn) -> None:
        """
        Add a reward function to the network.

        Args:
            name: Reward name
            reward_fn: Callable that takes RewardState and returns value
        """
        self._rewards[name] = reward_fn
        self._reset_struct()

    def get_reward(self, name: str):
        """Get reward function by name."""
        return self._rewards.get(name)

    def get_rewards(self) -> Dict[str, object]:
        """Get all rewards."""
        return self._rewards.copy()

    # =====================================================================
    # STATE INITIALIZATION METHODS
    # =====================================================================

    def init_default(self) -> None:
        """Initialize network with default state (all nodes empty)."""
        for node in self._nodes:
            if node.is_stateful():
                node.state = np.zeros(len(self._classes))

    def has_init_state(self) -> bool:
        """Check if network has initialized state."""
        return any(node.state is not None for node in self._nodes if node.is_stateful())

    def get_state(self) -> Dict[str, np.ndarray]:
        """
        Get current state of all stateful nodes.

        Returns:
            Dict mapping node names to state vectors
        """
        state = {}
        for node in self._nodes:
            if node.is_stateful() and node.state is not None:
                state[node.name] = node.state.copy()
        return state

    def set_state(self, state: Dict[str, np.ndarray]) -> None:
        """
        Set state for nodes.

        Args:
            state: Dict mapping node names to state vectors
        """
        for node_name, state_vec in state.items():
            node = self.get_node_by_name(node_name)
            if node is None:
                raise ValueError(f"[{self.name}] Node '{node_name}' not found")
            if not node.is_stateful():
                raise ValueError(f"[{self.name}] Node '{node_name}' is not stateful")
            node.state = state_vec

        self._reset_struct()

    def state(self) -> Dict[str, np.ndarray]:
        """
        Get current state of all stateful nodes (alias for get_state).

        Returns:
            Dict mapping node names to state vectors
        """
        return self.get_state()

    def init_from_marginal_and_started(self, n, s) -> None:
        """
        Initialize network state from marginal queue lengths and started jobs.

        Args:
            n: Marginal queue lengths matrix (nstations x nclasses).
               n[i][r] is the number of jobs of class r at station i.
            s: Started jobs matrix (nstations x nclasses).
               s[i][r] is the number of jobs of class r in service at station i.
        """
        n = np.atleast_2d(n)
        s = np.atleast_2d(s)

        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Validate dimensions
        if n.shape[0] != nstations:
            raise ValueError(f"n must have {nstations} rows (one per station)")
        if s.shape[0] != nstations:
            raise ValueError(f"s must have {nstations} rows (one per station)")

        # Set state for each stateful node
        for i, station in enumerate(self._stations):
            if hasattr(station, 'set_state'):
                # Combine queue length and started jobs into state
                # State format depends on scheduling strategy
                state_vec = np.zeros(nclasses)
                for k in range(min(nclasses, n.shape[1])):
                    state_vec[k] = n[i, k]
                station.set_state(state_vec)

        self._has_state = True
        self._reset_struct()

    # PascalCase alias
    initFromMarginalAndStarted = init_from_marginal_and_started

    # =====================================================================
    # UTILITY METHODS
    # =====================================================================

    def to_java(self):
        """Convert Network for JPype interop (future)."""
        raise NotImplementedError("Java interop not yet implemented")

    def __repr__(self) -> str:
        return (
            f"Network('{self.name}', "
            f"nodes={len(self._nodes)}, "
            f"stations={len(self._stations)}, "
            f"classes={len(self._classes)})"
        )

    # =====================================================================
    # COPY METHOD
    # =====================================================================

    def copy(self) -> 'Network':
        """
        Create a deep copy of the network.

        This method creates a new Network instance with copies of all nodes,
        classes, and routing. The copied network is independent of the original
        and can be modified without affecting the original.

        This is required for model transformations like Heidelberger-Trivedi (H-T)
        for fork-join networks.

        Returns:
            Network: A deep copy of this network

        Example:
            >>> model = Network('Original')
            >>> # ... build model ...
            >>> model_copy = model.copy()
            >>> # Modify model_copy without affecting model
        """
        import copy as copy_module

        # Create new network with same name
        new_network = Network(self.name)

        # Copy allow_replace flag if set
        if hasattr(self, 'allow_replace'):
            new_network.allow_replace = self.allow_replace

        # Build mapping from old nodes/classes to new ones
        node_map = {}  # old_node -> new_node
        class_map = {}  # old_class -> new_class

        # First pass: create copies of all nodes
        for old_node in self._nodes:
            # Deep copy the node
            new_node = copy_module.deepcopy(old_node)

            # Update model reference
            new_node._model = new_network

            # Reset index (will be set when added to network)
            new_node._node_index = len(new_network._nodes)

            # Add to new network's node list directly (bypass add_node to avoid re-validation)
            new_network._nodes.append(new_node)
            node_map[old_node] = new_node

            # Track stations separately
            if new_node.is_station():
                new_node._station_index = len(new_network._stations)
                new_network._stations.append(new_node)

        # Second pass: create copies of all classes
        for old_class in self._classes:
            # Deep copy the class
            new_class = copy_module.deepcopy(old_class)

            # Update index
            new_class._index = len(new_network._classes)

            # Update reference station to point to new node
            if hasattr(new_class, '_refstat') and new_class._refstat is not None:
                if new_class._refstat in node_map:
                    new_class._refstat = node_map[new_class._refstat]
            if hasattr(new_class, '_reference_station') and new_class._reference_station is not None:
                if new_class._reference_station in node_map:
                    new_class._reference_station = node_map[new_class._reference_station]

            new_network._classes.append(new_class)
            class_map[old_class] = new_class

        # Third pass: update node references (e.g., Fork references in Join nodes)
        for old_node, new_node in node_map.items():
            # Update Fork reference in Join nodes
            if hasattr(new_node, '_fork') and new_node._fork is not None:
                if new_node._fork in node_map:
                    new_node._fork = node_map[new_node._fork]

            # Update service/arrival process class references
            if hasattr(new_node, '_service_process') and new_node._service_process:
                new_service = {}
                for old_cls, dist in new_node._service_process.items():
                    if old_cls in class_map:
                        new_service[class_map[old_cls]] = dist
                    else:
                        new_service[old_cls] = dist
                new_node._service_process = new_service

            if hasattr(new_node, '_arrival_process') and new_node._arrival_process:
                new_arrival = {}
                for old_cls, dist in new_node._arrival_process.items():
                    if old_cls in class_map:
                        new_arrival[class_map[old_cls]] = dist
                    else:
                        new_arrival[old_cls] = dist
                new_node._arrival_process = new_arrival

            # Update routing strategies class references
            if hasattr(new_node, '_routing_strategies') and new_node._routing_strategies:
                new_strategies = {}
                for old_cls, strategy in new_node._routing_strategies.items():
                    if old_cls in class_map:
                        new_strategies[class_map[old_cls]] = strategy
                    else:
                        new_strategies[old_cls] = strategy
                new_node._routing_strategies = new_strategies

        # Copy routing matrix with updated references
        if self._routing_matrix is not None:
            new_rm = RoutingMatrix(new_network)
            old_routes = self._routing_matrix._routes
            for (old_from_class, old_to_class), routes in old_routes.items():
                new_from_class = class_map.get(old_from_class, old_from_class)
                new_to_class = class_map.get(old_to_class, old_to_class)
                for (old_from_node, old_to_node), prob in routes.items():
                    new_from_node = node_map.get(old_from_node, old_from_node)
                    new_to_node = node_map.get(old_to_node, old_to_node)
                    new_rm.set(new_from_class, new_to_class, new_from_node, new_to_node, prob)
            new_network._routing_matrix = new_rm

        # Copy rewards
        new_network._rewards = copy_module.deepcopy(self._rewards)

        # Reset cached values (they will be recomputed when needed)
        new_network._connections = None
        new_network._sn = None
        new_network._has_struct = False
        new_network._source_idx = -1
        new_network._sink_idx = -1

        return new_network

    # =====================================================================
    # CAMELCASE ALIASES FOR MATLAB/JAVA API COMPATIBILITY
    # =====================================================================

    # Link/connection methods
    addLink = add_link
    addLinks = add_links

    # Routing methods
    initRoutingMatrix = init_routing_matrix
    getRoutingMatrix = get_routing_matrix

    # Node/class accessors
    getNodes = get_nodes
    getStations = get_stations
    getClasses = get_classes
    getNumberOfNodes = get_number_of_nodes
    getNumberOfStations = get_number_of_stations
    getNumberOfClasses = get_number_of_classes
    getNodeByName = get_node_by_name
    getSource = get_source
    getSink = get_sink

    # State management
    refreshStruct = refresh_struct
    resetStruct = reset_struct
    getStruct = get_struct

    # Reward methods
    setReward = set_reward
    getReward = get_reward
    getRewards = get_rewards

    # Static routing helper
    @staticmethod
    def serial_routing(*args):
        """
        Create a routing probability matrix for serial (tandem) routing through nodes.

        Jobs flow from each node to the next in the provided order.
        For closed networks (last node is not a Sink), the last node
        automatically routes back to the first node to form a cycle.

        If a node appears multiple times in the sequence (e.g., for cyclic routing
        where the first node is repeated at the end), the duplicate is mapped back
        to the original node index to create a properly sized routing matrix.

        Args:
            *args: Either a single list of nodes, or nodes passed as separate arguments

        Returns:
            2D list of routing probabilities, where result[i][j] is the
            probability of routing from node i to node j.
        """
        # Handle both calling conventions:
        # serial_routing([node1, node2, node3]) - list form
        # serial_routing(node1, node2, node3) - variadic form
        # serial_routing(np.array([node1, node2, node3])) - numpy array form
        import numpy as np
        if len(args) == 1 and isinstance(args[0], (list, tuple, np.ndarray)):
            nodes = list(args[0])
        else:
            nodes = list(args)

        if len(nodes) == 0:
            raise ValueError("serial_routing requires at least one node")

        # Build mapping from node to unique index, preserving order of first occurrence
        unique_nodes = []
        node_to_idx = {}
        for node in nodes:
            if node not in node_to_idx:
                node_to_idx[node] = len(unique_nodes)
                unique_nodes.append(node)

        # Map each position in the sequence to its unique node index
        position_to_idx = [node_to_idx[node] for node in nodes]

        n = len(unique_nodes)
        # Create 2D list of routing probabilities
        P = [[0.0] * n for _ in range(n)]

        # Forward routing: node[i] -> node[i+1], mapped to unique indices
        for i in range(len(nodes) - 1):
            from_idx = position_to_idx[i]
            to_idx = position_to_idx[i + 1]
            P[from_idx][to_idx] = 1.0

        # Close the loop for non-sink networks (closed networks)
        # Only if the last node in the sequence is different from the first
        # (If they're the same, the loop is already closed by the above logic)
        from .nodes import Sink
        if len(nodes) > 0 and not isinstance(nodes[-1], Sink):
            last_idx = position_to_idx[-1]
            first_idx = position_to_idx[0]
            if last_idx != first_idx:
                P[last_idx][first_idx] = 1.0

        return P

    # PascalCase alias for MATLAB compatibility
    serialRouting = serial_routing

    # =====================================================================
    # FACTORY METHODS FOR CREATING STANDARD NETWORK TOPOLOGIES
    # =====================================================================

    @staticmethod
    def cyclic(N, D, strategy, S=None):
        """
        Create a cyclic queueing network with specified scheduling strategies.

        Creates a closed queueing network where jobs cycle through stations
        in a round-robin fashion (1 -> 2 -> ... -> M -> 1).

        Args:
            N: Population vector [1 x R] or list - number of jobs per class
            D: Service demand matrix [M x R] - service demands at each station per class
            strategy: List of scheduling strategies for each station [M]
                     (e.g., SchedStrategy.FCFS, SchedStrategy.PS, SchedStrategy.INF)
            S: Number of servers per station [M x 1] or list (default: 1 for each station)

        Returns:
            Network: Configured closed queueing network

        Example:
            >>> N = [10]  # 10 jobs of class 1
            >>> D = [[0.5], [1.0]]  # Service demands at 2 stations
            >>> strategy = [SchedStrategy.PS, SchedStrategy.FCFS]
            >>> model = Network.cyclic(N, D, strategy)

        References:
            MATLAB: matlab/src/lang/JNetwork.m
            Java: jar/src/main/kotlin/jline/lang/Network.java
        """
        from .nodes import Queue, Delay
        from .classes import ClosedClass
        from ..distributions import Exp

        # Convert to numpy arrays
        N = np.atleast_2d(N)
        D = np.atleast_2d(D)

        # Ensure N is a row vector [1 x R]
        if N.shape[0] > 1 and N.shape[1] == 1:
            N = N.T

        M = D.shape[0]  # Number of stations
        R = D.shape[1]  # Number of classes

        # Default: 1 server per station
        if S is None:
            S = np.ones((M, 1))
        else:
            S = np.atleast_2d(S)
            # Ensure S is column vector [M x 1]
            if S.shape[0] == 1 and S.shape[1] == M:
                S = S.T

        # Create the network
        model = Network("Model")
        nodes = []

        nqueues = 0
        ndelays = 0

        # Create stations based on scheduling strategy
        for i in range(M):
            if strategy[i] == SchedStrategy.INF:
                ndelays += 1
                node = Delay(model, f"Delay{ndelays}")
            else:
                nqueues += 1
                node = Queue(model, f"Queue{nqueues}", strategy[i])
                node.setNumberOfServers(int(S[i, 0]))
            nodes.append(node)

        # Create job classes (closed classes)
        jobclasses = []
        for r in range(R):
            # Reference station is the first node
            newclass = ClosedClass(model, f"Class{r + 1}", int(N[0, r]), nodes[0])
            jobclasses.append(newclass)

        # Set service processes
        for i in range(M):
            for r in range(R):
                demand = D[i, r]
                if demand > 0:
                    # Use Exp.fitMean equivalent: rate = 1/mean
                    nodes[i].setService(jobclasses[r], Exp(1.0 / demand))

        # Create routing matrix (cyclic: 1 -> 2 -> ... -> M -> 1)
        P = model.init_routing_matrix()
        for r in range(R):
            # Circulant routing for each class
            for i in range(M):
                next_i = (i + 1) % M
                P.set(jobclasses[r], jobclasses[r], nodes[i], nodes[next_i], 1.0)

        model.link(P)
        return model

    @staticmethod
    def cyclicPsInf(N, D, Z, S=None):
        """
        Create a cyclic network with Delay (INF) stations followed by PS queues.

        This creates a closed queueing network where:
        - The first MZ stations are Delays (infinite server, think time stations)
        - The remaining M stations are PS (processor sharing) queues

        Args:
            N: Population vector [1 x R] or list - number of jobs per class
            D: Service demand matrix [M x R] - demands at queue stations per class
            Z: Think time matrix [MZ x R] - think times at delay stations per class
               If all zeros, no delay stations are created.
            S: Number of servers per queue station [M x 1] or list (default: 1 each)

        Returns:
            Network: Configured closed queueing network with delays and PS queues

        Example:
            >>> N = [10]  # 10 jobs
            >>> D = [[1.0], [0.5]]  # Demands at 2 queues
            >>> Z = [[5.0]]  # 5 seconds think time at 1 delay
            >>> model = Network.cyclicPsInf(N, D, Z)

        References:
            MATLAB: matlab/src/lang/JNetwork.m
            Java: jar/src/main/kotlin/jline/lang/Network.java
        """
        # Convert to numpy arrays
        D = np.atleast_2d(D)
        Z = np.atleast_2d(Z)

        M = D.shape[0]  # Number of queue stations
        MZ = Z.shape[0]  # Number of delay stations
        R = D.shape[1]  # Number of classes

        # If Z is all zeros, don't create delay stations
        if np.max(Z) == 0:
            MZ = 0

        # Build scheduling strategy array: INF for delays, PS for queues
        strategy = []
        for i in range(MZ):
            strategy.append(SchedStrategy.INF)
        for i in range(M):
            strategy.append(SchedStrategy.PS)

        # Build combined demand matrix [Z; D]
        if MZ > 0:
            Dnew = np.vstack([Z, D])
        else:
            Dnew = D

        # Build combined server count [inf for delays; S for queues]
        if S is None:
            S = np.ones((M, 1))
        else:
            S = np.atleast_2d(S)
            if S.shape[0] == 1 and S.shape[1] == M:
                S = S.T

        if MZ > 0:
            Sinf = np.full((MZ, 1), np.inf)
            Snew = np.vstack([Sinf, S])
        else:
            Snew = S

        return Network.cyclic(N, Dnew, strategy, Snew)

    @staticmethod
    def cyclicFcfs(N, D, S=None):
        """
        Create a cyclic queueing network with FCFS scheduling at all stations.

        Args:
            N: Population vector [1 x R] or list - number of jobs per class
            D: Service demand matrix [M x R] - service demands at each station per class
            S: Number of servers per station [M x 1] or list (default: 1 each)

        Returns:
            Network: Configured closed queueing network with FCFS scheduling

        Example:
            >>> N = [10]  # 10 jobs
            >>> D = [[0.5], [1.0]]  # Service demands at 2 stations
            >>> model = Network.cyclicFcfs(N, D)

        References:
            MATLAB: matlab/src/lang/JNetwork.m
            Java: jar/src/main/kotlin/jline/lang/Network.java
        """
        D = np.atleast_2d(D)
        M = D.shape[0]

        # All FCFS scheduling
        strategy = [SchedStrategy.FCFS] * M

        return Network.cyclic(N, D, strategy, S)

    @staticmethod
    def cyclicPs(N, D, S=None):
        """
        Create a cyclic queueing network with PS scheduling at all stations.

        Args:
            N: Population vector [1 x R] or list - number of jobs per class
            D: Service demand matrix [M x R] - service demands at each station per class
            S: Number of servers per station [M x 1] or list (default: 1 each)

        Returns:
            Network: Configured closed queueing network with PS scheduling

        Example:
            >>> N = [10]  # 10 jobs
            >>> D = [[0.5], [1.0]]  # Service demands at 2 stations
            >>> model = Network.cyclicPs(N, D)

        References:
            MATLAB: matlab/src/lang/JNetwork.m
            Java: jar/src/main/kotlin/jline/lang/Network.java
        """
        D = np.atleast_2d(D)
        M = D.shape[0]

        # All PS scheduling
        strategy = [SchedStrategy.PS] * M

        return Network.cyclic(N, D, strategy, S)

    # =====================================================================
    # VISUALIZATION METHODS
    # =====================================================================

    def jsimgView(self) -> bool:
        """
        Open the model in JMT's JSIMgraph graphical editor.

        This method exports the network to JSIMG format (JMT simulation model)
        and opens it in JSIMgraph for viewing and editing.

        Returns:
            True if JMT was launched successfully, False otherwise

        Raises:
            ImportError: If io module is not available

        Example:
            >>> model = Network('MyModel')
            >>> # ... build model ...
            >>> model.jsimgView()  # Opens in JMT graphical editor

        References:
            MATLAB: matlab/src/lang/@MNetwork/jsimgView.m
        """
        import tempfile
        import os
        from ..api.io import jsimg_view, line_printf
        from ..api.solvers.jmt.handler import _write_jsim_file, SolverJMTOptions

        # Compile model if needed
        if not self._has_struct:
            self.link(self._routing_matrix)

        # Get NetworkStruct
        sn = self.get_struct()

        # Create temp file for JSIMG
        fd, jsimg_file = tempfile.mkstemp(suffix='.jsimg', prefix='model_')
        os.close(fd)

        # Write model to JSIMG format using the proper handler
        options = SolverJMTOptions()
        _write_jsim_file(sn, jsimg_file, options)

        line_printf('JMT Model: %s\n', jsimg_file)

        # Open in JSIMgraph
        return jsimg_view(jsimg_file)

    # snake_case alias
    jsimg_view = jsimgView

    def jsimwView(self) -> bool:
        """
        Open the model in JMT's JSIMwiz wizard interface.

        This method exports the network to JSIMG format and opens it
        in JSIMwiz for wizard-style configuration and simulation.

        Returns:
            True if JMT was launched successfully, False otherwise

        Example:
            >>> model = Network('MyModel')
            >>> # ... build model ...
            >>> model.jsimwView()  # Opens in JMT wizard

        References:
            MATLAB: matlab/src/lang/@MNetwork/jsimwView.m
        """
        import tempfile
        import os
        from ..api.io import jsimw_view, line_printf
        from ..api.solvers.jmt.handler import _write_jsim_file, SolverJMTOptions

        # Compile model if needed
        if not self._has_struct:
            self.link(self._routing_matrix)

        # Get NetworkStruct
        sn = self.get_struct()

        # Create temp file for JSIMG
        fd, jsimg_file = tempfile.mkstemp(suffix='.jsimg', prefix='model_')
        os.close(fd)

        # Write model to JSIMG format using the proper handler
        options = SolverJMTOptions()
        _write_jsim_file(sn, jsimg_file, options)

        line_printf('JMT Model: %s\n', jsimg_file)

        # Open in JSIMwiz
        return jsimw_view(jsimg_file)

    # snake_case alias
    jsimw_view = jsimwView

    def view(self) -> bool:
        """
        Open the model in JMT's graphical editor (alias for jsimgView).

        This is a convenience alias that opens the model in JSIMgraph,
        providing a visual representation of the queueing network.

        Returns:
            True if JMT was launched successfully, False otherwise

        Example:
            >>> model = Network('MyModel')
            >>> # ... build model ...
            >>> model.view()  # Opens in JMT graphical editor

        References:
            MATLAB: matlab/src/lang/@MNetwork/view.m
        """
        return self.jsimgView()

    def modelView(self) -> None:
        """
        Open the model in ModelVisualizer.

        Note: ModelVisualizer is a Java Swing GUI component that is not
        available in the Python native implementation. Use jsimgView() or
        view() for visualization, which opens the model in JMT's graphical
        editor.

        Raises:
            NotImplementedError: ModelVisualizer requires Java GUI

        References:
            MATLAB: matlab/src/lang/@MNetwork/modelView.m
        """
        raise NotImplementedError(
            "modelView() requires Java GUI (ModelVisualizer). "
            "Use view() or jsimgView() to open the model in JMT's graphical editor."
        )

    # snake_case alias
    model_view = modelView

    # =====================================================================
    # STATIC FACTORY METHODS
    # =====================================================================

    @staticmethod
    def tandemPsInf(lambda_rates: np.ndarray, D: np.ndarray,
                    Z: Optional[np.ndarray] = None) -> 'Network':
        """
        Create a tandem network with PS queues and INF (delay) stations.

        Creates an open queueing network with:
        - Source node generating arrivals
        - Optional Delay stations (INF scheduling) from Z
        - Queue stations (PS scheduling) from D
        - Sink node

        Args:
            lambda_rates: Array of arrival rates for each class (R,)
            D: Service time matrix for queues (M x R), where M is number
               of PS queues and R is number of classes
            Z: Optional delay service times (Mz x R), where Mz is number
               of delay stations

        Returns:
            Network: Configured tandem queueing network

        Example:
            >>> lambda_rates = np.array([0.02, 0.04])
            >>> D = np.array([[10, 5], [5, 9]])
            >>> Z = np.array([[91, 92]])
            >>> model = Network.tandemPsInf(lambda_rates, D, Z)

        References:
            MATLAB: matlab/src/lang/@MNetwork/tandemPsInf.m
        """
        from .nodes import Source, Queue, Delay, Sink
        from .classes import OpenClass
        from ..distributions import Exp

        if Z is None:
            Z = np.array([]).reshape(0, D.shape[1]) if len(D.shape) > 1 else np.array([])

        # Ensure Z and D are 2D
        Z = np.atleast_2d(Z) if Z.size > 0 else np.zeros((0, D.shape[1]))
        D = np.atleast_2d(D)

        M = D.shape[0]  # Number of PS queues
        Mz = Z.shape[0]  # Number of delay stations
        R = D.shape[1]  # Number of classes

        # Build strategy list: INF for delays, PS for queues
        strategies = [SchedStrategy.INF] * Mz + [SchedStrategy.PS] * M

        # Combine service times
        if Mz > 0:
            combined_D = np.vstack([Z, D])
        else:
            combined_D = D

        return Network.tandem(lambda_rates, combined_D, strategies)

    @staticmethod
    def tandem(lambda_rates: np.ndarray, D: np.ndarray,
               strategies: List) -> 'Network':
        """
        Create a tandem network with specified scheduling strategies.

        Creates an open queueing network in tandem configuration.

        Args:
            lambda_rates: Array of arrival rates for each class (R,)
            D: Service time matrix (M x R), where D[i,r] is mean service
               time of class r at station i
            strategies: List of scheduling strategies for each station

        Returns:
            Network: Configured tandem queueing network

        References:
            MATLAB: matlab/src/lang/@MNetwork/tandem.m
        """
        from .nodes import Source, Queue, Delay, Sink
        from .classes import OpenClass
        from ..distributions import Exp

        D = np.atleast_2d(D)
        M, R = D.shape

        model = Network('Model')

        # Create nodes
        nodes = []
        source = Source(model, 'Source')
        nodes.append(source)

        for i in range(M):
            strategy = strategies[i] if i < len(strategies) else SchedStrategy.FCFS
            if strategy == SchedStrategy.INF:
                node = Delay(model, f'Delay{i+1}')
            else:
                node = Queue(model, f'Queue{i+1}', strategy)
            nodes.append(node)

        sink = Sink(model, 'Sink')
        nodes.append(sink)

        # Create job classes
        jobclasses = []
        for r in range(R):
            jobclass = OpenClass(model, f'Class{r+1}')
            jobclasses.append(jobclass)

        # Set arrival and service rates
        lambda_arr = np.atleast_1d(lambda_rates)
        for r in range(R):
            # Arrival rate at source
            arr_rate = lambda_arr[r] if r < len(lambda_arr) else 1.0
            source.set_arrival(jobclasses[r], Exp.fit_mean(1.0 / arr_rate))

            # Service times at each station
            for i in range(M):
                service_time = D[i, r] if D[i, r] > 0 else 1.0
                nodes[i + 1].set_service(jobclasses[r], Exp.fit_mean(service_time))

        # Create serial routing
        P = model.init_routing_matrix()
        for r in range(R):
            for i in range(len(nodes) - 1):
                P.set(jobclasses[r], jobclasses[r], nodes[i], nodes[i + 1], 1.0)
        model.link(P)

        return model

    # Snake_case aliases
    tandem_ps_inf = tandemPsInf
    cyclic_ps_inf = cyclicPsInf
    cyclic_fcfs = cyclicFcfs
    cyclic_ps = cyclicPs


__all__ = ['Network']
