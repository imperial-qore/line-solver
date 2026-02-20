"""
Network class for LINE native Python implementation.

This module provides the Network class for constructing and managing queueing networks.
Ported from MATLAB implementation in matlab/src/lang/@MNetwork/
"""

from typing import Optional, Dict, List, Tuple, Union
import numpy as np

from .base import (
    Element, ElementType, Network as NetworkBase, NetworkElement, Node, StatefulNode, Station,
    JobClass, JobClassType, NodeType, SchedStrategy, RoutingStrategy, DropStrategy
)
from .region import Region
from .routing import RoutingMatrix
from ..api.sn.network_struct import NetworkStruct, DropStrategy as SnDropStrategy
from ..api.sn.utils import sn_print_routing_matrix
from ..constants import ProcessType, GlobalConstants


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
        self._links = set()  # Set of (source_idx, dest_idx) tuples for explicit links
        self._connections = None  # Adjacency matrix (lazy-computed)
        self._routing_matrix = None  # Routing matrix
        self._sn = None  # Compiled NetworkStruct
        self._has_struct = False  # Whether struct is valid
        self._rewards = {}  # Dict[reward_name, reward_function]
        self._source_idx = -1  # Cached source node index
        self._sink_idx = -1  # Cached sink node index
        self._log_path = None  # Path for logger output files
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

        # Note: Unlike the old implementation that added probability 1.0 for all classes,
        # we now only track connectivity here. Routing probabilities are computed later
        # in _refresh_routing() based on the routing strategy (PROB vs RAND).
        # This matches MATLAB's addLink behavior which only sets the connection matrix.

        # Track the link explicitly
        src_idx = source._node_index if hasattr(source, '_node_index') else self._nodes.index(source)
        dst_idx = dest._node_index if hasattr(dest, '_node_index') else self._nodes.index(dest)
        self._links.add((src_idx, dst_idx))

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
            # Handle nested list format P[r][s], P[r], or simple 2D list
            rm = self.init_routing_matrix()
            nodes = self.get_nodes()
            classes = self.get_classes()
            nclasses = len(classes)
            nnodes = len(nodes)

            # Determine the format of the routing matrix:
            # 1. P[r][s] format: routing_matrix[r][s] is a 2D node routing matrix (class switching)
            # 2. P[r] format: routing_matrix[r] is a 2D node routing matrix (per-class routing, no switching)
            # 3. Simple 2D list: same routing for all classes

            routing_format = 'simple'  # default

            if len(routing_matrix) == nclasses and nclasses > 0:
                first_elem = routing_matrix[0]
                if isinstance(first_elem, (list, np.ndarray)):
                    # Check if first_elem is a list of length nclasses with 2D matrices
                    # This is P[r][s] format (class switching)
                    if isinstance(first_elem, list) and len(first_elem) == nclasses:
                        # Check if any first_elem[x] is a 2D matrix (list of lists or 2D array)
                        # Find first non-None entry to test format
                        test_entry = None
                        for entry in first_elem:
                            if entry is not None:
                                test_entry = entry
                                break
                        if test_entry is None:
                            # Try other rows for a non-None entry
                            for r_check in range(nclasses):
                                for s_check in range(nclasses):
                                    if routing_matrix[r_check][s_check] is not None:
                                        test_entry = routing_matrix[r_check][s_check]
                                        break
                                if test_entry is not None:
                                    break
                        if test_entry is not None and isinstance(test_entry, (list, np.ndarray)):
                            test_arr = np.asarray(test_entry)
                            if test_arr.ndim == 2:
                                routing_format = 'class_switching'

                    # If not class_switching, check for per_class format
                    if routing_format == 'simple':
                        first_arr = np.asarray(first_elem)
                        if first_arr.ndim == 2 and first_arr.shape[0] == nnodes and first_arr.shape[1] == nnodes:
                            # P[r] format - each P[r] is a 2D node routing matrix
                            routing_format = 'per_class'

            if routing_format == 'class_switching':
                # P[r][s] format - routing from class r to class s
                for r in range(nclasses):
                    for s in range(nclasses):
                        if routing_matrix[r][s] is None:
                            continue
                        matrix_arr = np.asarray(routing_matrix[r][s])
                        if matrix_arr.size == 0:
                            continue
                        for i in range(min(nnodes, matrix_arr.shape[0])):
                            for j in range(min(nnodes, matrix_arr.shape[1])):
                                if matrix_arr[i, j] > 0:
                                    rm.set(classes[r], classes[s], nodes[i], nodes[j], matrix_arr[i, j])
            elif routing_format == 'per_class':
                # P[r] format - each class has its own routing matrix (no class switching)
                for r in range(nclasses):
                    matrix_arr = np.asarray(routing_matrix[r])
                    for i in range(min(nnodes, matrix_arr.shape[0])):
                        for j in range(min(nnodes, matrix_arr.shape[1])):
                            if matrix_arr[i, j] > 0:
                                rm.set(classes[r], classes[r], nodes[i], nodes[j], matrix_arr[i, j])
            else:
                # Simple 2D list - apply same routing to all classes
                matrix_arr = np.asarray(routing_matrix)
                for jobclass in classes:
                    for i in range(min(nnodes, matrix_arr.shape[0])):
                        for j in range(min(nnodes, matrix_arr.shape[1])):
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

        # Sync routing probabilities to nodes' _prob_routing attribute
        # This ensures routing strategy is set to PROB when explicit probabilities are used
        # (matches MATLAB's setProbRouting behavior in link.m)
        self._sync_prob_routing_from_matrix()

        # Insert auto-generated ClassSwitch nodes for class switching routes
        self._insert_auto_class_switches()

        self._refresh_routing_matrix()
        self._reset_struct()

        # Check for reducible routing (absorbing states)
        # This matches MATLAB's link.m lines 304-312
        # Use the RoutingMatrix directly to avoid triggering refresh_struct()
        # which can cause infinite recursion for models with Fork/Join nodes.
        if self._routing_matrix is not None:
            try:
                rt = self._routing_matrix.toMatrix()
                nstations = self.get_number_of_stations()
                nclasses = len(self._classes)
                if rt.size > 0 and nstations > 0 and nclasses > 0:
                    # Build adjacency from station-class routing matrix
                    adj = np.zeros((nstations, nstations))
                    for i in range(nstations):
                        for j in range(nstations):
                            block = rt[i*nclasses:(i+1)*nclasses, j*nclasses:(j+1)*nclasses]
                            if np.any(block > 0):
                                adj[i, j] = 1.0
                    # Quick check: any station with no outgoing edges (excluding self)?
                    has_absorbing = False
                    for i in range(nstations):
                        outgoing = adj[i, :].copy()
                        outgoing[i] = 0
                        if np.sum(outgoing) < 1e-10 and np.sum(adj[i, :]) < 1e-10:
                            has_absorbing = True
                            break
                    if has_absorbing:
                        import warnings
                        warnings.warn(
                            "[link] Reducible network topology detected, results may be unreliable.",
                            UserWarning
                        )
            except Exception:
                pass  # Skip ergodic check if routing matrix conversion fails

    def is_routing_ergodic(self, P=None):
        """
        Check if the queueing network routing matrix is ergodic (irreducible).

        This checks only the routing structure, not the full CTMC state space.
        A routing is ergodic if all stations communicate, meaning the routing
        matrix does not create absorbing states or disconnected components.

        Args:
            P: Optional routing matrix (list of lists). If not provided, it will
               be retrieved via get_linked_routing_matrix().

        Returns:
            Tuple of (is_ergodic, info) where info is a dict containing:
            - absorbingStations: list of station names that are absorbing
            - transientStations: list of station names that are transient
            - numSCCs: number of strongly connected components
            - isReducible: True if the routing creates a reducible structure
        """
        from scipy.sparse.csgraph import connected_components
        from scipy.sparse import csr_matrix

        info = {
            'absorbingStations': [],
            'transientStations': [],
            'numSCCs': 1,
            'isReducible': False
        }

        # Get routing matrix if not provided
        if P is None:
            P = self.get_linked_routing_matrix()
        if P is None or len(P) == 0:
            return True, info

        nclasses = len(self._classes)
        nnodes = len(self._nodes)
        node_names = [n.name for n in self._nodes]

        # Build aggregate adjacency matrix across all classes
        adj_matrix = np.zeros((nnodes, nnodes))
        for r in range(nclasses):
            for s in range(nclasses):
                if P[r][s] is not None:
                    p_rs = np.asarray(P[r][s])
                    if p_rs.size > 0:
                        rows = min(nnodes, p_rs.shape[0])
                        cols = min(nnodes, p_rs.shape[1])
                        adj_matrix[:rows, :cols] += (p_rs[:rows, :cols] > 0).astype(float)
        adj_matrix = (adj_matrix > 0).astype(float)

        # Find strongly connected components
        n_components, labels = connected_components(
            csr_matrix(adj_matrix), directed=True, connection='strong'
        )
        info['numSCCs'] = n_components

        # Check for absorbing stations (stations with only self-loops or no outgoing edges)
        station_idxs = [i for i, n in enumerate(self._nodes) if n in self._stations]

        for idx in station_idxs:
            # Check if this station has any outgoing transitions to other nodes
            outgoing = adj_matrix[idx, :].copy()
            outgoing[idx] = 0  # Exclude self-loop

            if np.sum(outgoing) == 0:
                # No outgoing transitions except possibly self-loop
                # Check if there's routing defined for this node
                has_routing = False
                for r in range(nclasses):
                    for s in range(nclasses):
                        if P[r][s] is not None:
                            p_rs = np.asarray(P[r][s])
                            if idx < p_rs.shape[0] and np.any(p_rs[idx, :] > 0):
                                has_routing = True
                                break
                    if has_routing:
                        break

                if has_routing:
                    info['absorbingStations'].append(node_names[idx])

        # Identify transient stations (stations in transient SCCs)
        # A SCC is recurrent if it has outgoing edges only to itself
        for i in range(n_components):
            scc_nodes = np.where(labels == i)[0]
            # Check if this SCC has outgoing edges to other SCCs
            scc_outgoing = False
            for node in scc_nodes:
                for other_node in range(nnodes):
                    if labels[other_node] != i and adj_matrix[node, other_node] > 0:
                        scc_outgoing = True
                        break
                if scc_outgoing:
                    break

            if scc_outgoing:
                # This SCC is transient
                for node_idx in scc_nodes:
                    if node_idx in station_idxs:
                        node_name = node_names[node_idx]
                        if node_name not in info['absorbingStations']:
                            info['transientStations'].append(node_name)

        # Remove Sink from absorbing list (it's expected to be absorbing in open networks)
        sink = self.get_sink()
        if sink is not None:
            sink_name = sink.name
            info['absorbingStations'] = [s for s in info['absorbingStations'] if s != sink_name]

        # Determine ergodicity
        has_absorbing_stations = len(info['absorbingStations']) > 0
        # Check for multiple recurrent SCCs
        recurrent_scc_count = 0
        for i in range(n_components):
            scc_nodes = np.where(labels == i)[0]
            # Check if this SCC has outgoing edges to other SCCs
            scc_outgoing = False
            for node in scc_nodes:
                for other_node in range(nnodes):
                    if labels[other_node] != i and adj_matrix[node, other_node] > 0:
                        scc_outgoing = True
                        break
                if scc_outgoing:
                    break
            if not scc_outgoing:
                recurrent_scc_count += 1

        has_multiple_recurrent_sccs = recurrent_scc_count > 1

        if has_absorbing_stations or has_multiple_recurrent_sccs:
            is_ergodic = False
            info['isReducible'] = True
        else:
            is_ergodic = True
            info['isReducible'] = False

        return is_ergodic, info

    def _sync_prob_routing_from_matrix(self) -> None:
        """
        Sync routing probabilities from _routing_matrix to nodes' _prob_routing attribute.

        This ensures that when link() is called with explicit routing probabilities,
        each source node's _prob_routing is set so that refresh_struct() will
        correctly set the routing strategy to PROB instead of RAND.

        This matches MATLAB's behavior in link.m where setProbRouting is called
        for each non-zero routing probability (line 274).
        """
        if self._routing_matrix is None:
            return

        # Clear existing _prob_routing on all nodes to avoid stale routes
        # (e.g., when link() is called again after link_and_log() copies the model)
        for node in self.get_nodes():
            if hasattr(node, '_prob_routing'):
                node._prob_routing = {}

        # Iterate through all routes in the routing matrix
        for (class_src, class_dst), routes in self._routing_matrix._routes.items():
            for (node_src, node_dst), prob in routes.items():
                if prob > 0:
                    # Set _prob_routing on the source node
                    # This ensures routing strategy will be PROB instead of RAND
                    node_src.set_prob_routing(class_src, node_dst, prob)

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
            # Handle integer indices (from P[0] = ... syntax) or JobClass objects
            if isinstance(class_src, (int, np.integer)):
                src_class_idx = class_src
            else:
                src_class_idx = class_src._index if hasattr(class_src, '_index') else classes.index(class_src)
            if isinstance(class_dst, (int, np.integer)):
                dst_class_idx = class_dst
            else:
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
            # Normalize rows and add identity for classes with no switching
            # (matches MATLAB link.m lines 191-197)
            for r in range(nclasses):
                row_sum = np.sum(csmat[r, :])
                if row_sum > FINE_TOL:
                    csmat[r, :] = csmat[r, :] / row_sum
                else:
                    # Class r has no switching through this edge - add identity (pass-through)
                    csmat[r, r] = 1.0

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
            # Handle integer indices (from P[0] = ... syntax) or JobClass objects
            if isinstance(class_src, (int, np.integer)):
                src_class_idx = class_src
            else:
                src_class_idx = class_src._index if hasattr(class_src, '_index') else classes.index(class_src)
            if isinstance(class_dst, (int, np.integer)):
                dst_class_idx = class_dst
            else:
                dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else classes.index(class_dst)

            for (node_src, node_dst), prob in routes.items():
                src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else nodes.index(node_src)
                dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else nodes.index(node_dst)

                # Route through ClassSwitch if there's a CS node on this edge
                # ALL traffic through a CS edge should go through the CS node - the switch
                # matrix handles class transformation (including same-class pass-through)
                if prob > 0 and (src_node_idx, dst_node_idx) in cs_node_map:
                    # Route through CS node
                    cs_node = cs_node_map[(src_node_idx, dst_node_idx)]

                    # P[r,r](src, CS) += P[r,s](src, dst) - accumulate all traffic to CS
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
                    # Keep original route (no CS node on this edge)
                    key = (class_src, class_dst)
                    if key not in new_routes:
                        new_routes[key] = {}
                    route_key = (node_src, node_dst)
                    new_routes[key][route_key] = new_routes[key].get(route_key, 0) + prob

        # Save original routes before modifying (used by toMatrix() for rt computation)
        # Shallow copy the structure but preserve references to node/class objects
        original_routes = {}
        for key, routes_dict in self._routing_matrix._routes.items():
            original_routes[key] = dict(routes_dict)  # Shallow copy of inner dict
        self._routing_matrix._original_routes = original_routes

        # Update routing matrix with ClassSwitch node routes
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

    def get_log_path(self) -> str:
        """
        Get the path for logger output files.

        Returns:
            Path for log files, or temp directory if not set
        """
        import tempfile
        if self._log_path is None:
            return tempfile.gettempdir()
        return self._log_path

    def set_log_path(self, path: str) -> None:
        """
        Set the path for logger output files.

        Args:
            path: Directory path for log files
        """
        self._log_path = path

    # CamelCase aliases
    getLogPath = get_log_path
    setLogPath = set_log_path

    def link_and_log(self, P, is_node_logged: list, log_path: str = None):
        """
        Link the network with logging enabled for specified nodes.

        This method modifies the network by inserting Logger nodes before
        and after the logged nodes to capture arrival and departure timestamps.

        Ported from MATLAB's MNetwork.linkAndLog.

        Args:
            P: Routing matrix (dict of dicts of numpy arrays)
            is_node_logged: Boolean list indicating which nodes should be logged
            log_path: Path for log files (optional, uses model's log path if not set)

        Returns:
            Tuple of (loggerBefore, loggerAfter) - lists of Logger nodes created
        """
        import os
        from .nodes import Logger

        self.reset_struct()

        if self._has_struct:
            self.reset_network()

        if log_path is not None:
            self.set_log_path(log_path)
        else:
            log_path = self.get_log_path()

        R = self.get_number_of_classes()
        M_nodes = self.get_number_of_nodes()

        if len(is_node_logged) != M_nodes:
            raise ValueError(f"is_node_logged length ({len(is_node_logged)}) must match number of nodes ({M_nodes})")

        is_node_logged = list(is_node_logged)

        # Don't log Source or Sink
        source = self.get_source()
        if source is not None:
            sink_idx = self.get_index_sink_node()
            if sink_idx >= 0 and sink_idx < len(is_node_logged) and is_node_logged[sink_idx]:
                is_node_logged[sink_idx] = False

            source_idx = self.get_index_source_node()
            if source_idx >= 0 and source_idx < len(is_node_logged) and is_node_logged[source_idx]:
                is_node_logged[source_idx] = False

        logger_before = []
        logger_after = []

        # Create departure loggers FIRST (to match JAR node ordering)
        for ind in range(M_nodes):
            if is_node_logged[ind]:
                node_name = self._nodes[ind].name
                log_file = os.path.join(log_path, f"{node_name}-Dep.csv")
                logger = Logger(self, f"Dep_{node_name}", log_file)
                for r in range(R):
                    logger.set_routing(self._classes[r], 1)  # RAND routing
                logger_after.append(logger)

        # Create arrival loggers SECOND
        for ind in range(M_nodes):
            if is_node_logged[ind]:
                node_name = self._nodes[ind].name
                log_file = os.path.join(log_path, f"{node_name}-Arv.csv")
                logger = Logger(self, f"Arv_{node_name}", log_file)
                for r in range(R):
                    logger.set_routing(self._classes[r], 1)  # RAND routing
                logger_before.append(logger)

        # Build new routing matrix with loggers
        # Original nodes: 0..M_nodes-1
        # Logger after (Dep): M_nodes..M_nodes+len(logged)-1  (created first)
        # Logger before (Arv): M_nodes+len(logged)..  (created second)
        M_nodes_new = 3 * M_nodes  # Upper bound
        logged_indices = [i for i, logged in enumerate(is_node_logged) if logged]
        num_logged = len(logged_indices)

        # Build routing mapping
        # For each (r, s) class pair, create new routing matrix
        new_P = {}
        for r in range(R):
            class_r = self._classes[r]
            new_P[class_r] = {}
            for s in range(R):
                class_s = self._classes[s]
                new_routing = np.zeros((M_nodes + 2 * num_logged, M_nodes + 2 * num_logged))

                # Get original routing
                # P can be either:
                # 1. List of lists: P[r_idx][s_idx] is routing matrix
                # 2. Dict of dicts: P[class_r][class_s] is routing matrix
                if isinstance(P, list):
                    # List format - use integer indices
                    if r < len(P) and s < len(P[r]) and P[r][s] is not None:
                        orig_routing = P[r][s]
                    else:
                        orig_routing = np.zeros((M_nodes, M_nodes))
                elif isinstance(P, dict):
                    # Dict format - use class objects
                    if class_r in P and class_s in P[class_r]:
                        orig_routing = P[class_r][class_s]
                    else:
                        orig_routing = np.zeros((M_nodes, M_nodes))
                else:
                    orig_routing = np.zeros((M_nodes, M_nodes))

                # Transform routing based on logging
                # Dep loggers at M_nodes+k, Arv loggers at M_nodes+num_logged+k
                for ind in range(M_nodes):
                    for jnd in range(M_nodes):
                        if orig_routing[ind, jnd] > 0:
                            i_logged = is_node_logged[ind]
                            j_logged = is_node_logged[jnd]

                            if i_logged and j_logged:
                                # Link loggerDep_i to loggerArv_j
                                i_dep = M_nodes + logged_indices.index(ind)
                                j_arv = M_nodes + num_logged + logged_indices.index(jnd)
                                new_routing[i_dep, j_arv] = orig_routing[ind, jnd]
                            elif i_logged and not j_logged:
                                # Link loggerDep_i to j
                                i_dep = M_nodes + logged_indices.index(ind)
                                new_routing[i_dep, jnd] = orig_routing[ind, jnd]
                            elif not i_logged and j_logged:
                                # Link i to loggerArv_j
                                j_arv = M_nodes + num_logged + logged_indices.index(jnd)
                                new_routing[ind, j_arv] = orig_routing[ind, jnd]
                            else:
                                # Link i to j (no loggers)
                                new_routing[ind, jnd] = orig_routing[ind, jnd]

                # Add logger chain links (for same class only)
                if r == s:
                    for k, ind in enumerate(logged_indices):
                        dep_idx = M_nodes + k
                        arv_idx = M_nodes + num_logged + k
                        # loggerArv -> node
                        new_routing[arv_idx, ind] = 1.0
                        # node -> loggerDep
                        new_routing[ind, dep_idx] = 1.0

                new_P[class_r][class_s] = new_routing

        # Now reduce the routing matrix to only include used nodes
        used_nodes = list(range(M_nodes))  # Original nodes
        for k in range(num_logged):
            used_nodes.append(M_nodes + k)  # Dep loggers (created first)
        for k in range(num_logged):
            used_nodes.append(M_nodes + num_logged + k)  # Arv loggers (created second)

        # Convert nested dict format to flat tuple-key dict format expected by link()
        # {(from_class, to_class): matrix}
        final_P = {}
        for class_r in new_P:
            for class_s in new_P[class_r]:
                full_matrix = new_P[class_r][class_s]
                reduced = full_matrix[np.ix_(used_nodes, used_nodes)]
                final_P[(class_r, class_s)] = reduced

        self.link(final_P)
        return logger_before, logger_after

    # CamelCase alias
    linkAndLog = link_and_log

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
        """Compute adjacency matrix from explicit links and routing matrix."""
        nnodes = len(self._nodes)
        connections = np.zeros((nnodes, nnodes), dtype=int)

        # Add explicit links from add_link() calls
        for (src_idx, dst_idx) in self._links:
            if 0 <= src_idx < nnodes and 0 <= dst_idx < nnodes:
                connections[src_idx, dst_idx] = 1

        # Also add any routes from routing_matrix (e.g., from set_prob_routing)
        if self._routing_matrix is not None:
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                for (node_src, node_dst), prob in routes.items():
                    if prob > 0:
                        src_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                        dst_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)
                        if 0 <= src_idx < nnodes and 0 <= dst_idx < nnodes:
                            connections[src_idx, dst_idx] = 1

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

    @property
    def classes(self) -> List[JobClass]:
        """Get list of all job classes (property for compatibility with MATLAB API)."""
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
                # Get population from ClosedClass
                if hasattr(jobclass, 'getNumberOfJobs'):
                    njobs[i] = jobclass.getNumberOfJobs()
                elif hasattr(jobclass, '_njobs'):
                    njobs[i] = jobclass._njobs
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

    def get_index_source_node(self) -> int:
        """Get the index of the Source node (-1 if not present)."""
        if self._source_idx >= 0:
            return self._source_idx

        for i, node in enumerate(self._nodes):
            if node.node_type == NodeType.SOURCE:
                self._source_idx = i
                return i

        return -1

    def get_index_sink_node(self) -> int:
        """Get the index of the Sink node (-1 if not present)."""
        if self._sink_idx >= 0:
            return self._sink_idx

        for i, node in enumerate(self._nodes):
            if node.node_type == NodeType.SINK:
                self._sink_idx = i
                return i

        return -1

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

    def init_used_features(self):
        """Initialize the used features set."""
        from ..solvers.base import SolverFeatureSet
        self._used_features = SolverFeatureSet()

    def set_used_lang_feature(self, feature: str):
        """Set a feature as used in this model."""
        if not hasattr(self, '_used_features') or self._used_features is None:
            self.init_used_features()
        self._used_features.set_true(feature)

    def get_used_lang_features(self):
        """
        Get the features used by this model.

        Returns:
            SolverFeatureSet: The set of features used by this model.
        """
        from ..solvers.base import SolverFeatureSet
        from .base import SchedStrategy, RoutingStrategy, ReplacementStrategy, JobClassType

        self.init_used_features()

        # Check for open and closed classes
        for jobclass in self._classes:
            if jobclass.jobclass_type == JobClassType.OPEN:
                self.set_used_lang_feature('OpenClass')
            elif jobclass.jobclass_type == JobClassType.CLOSED:
                self.set_used_lang_feature('ClosedClass')

        # Check features from nodes
        from .nodes import Queue, Source, Sink, Cache, ClassSwitch, Fork, Join, Router, Place, Transition

        for node in self._nodes:
            node_type = type(node).__name__

            if node_type == 'Queue':
                self.set_used_lang_feature('Queue')
                # Check scheduling strategy
                if hasattr(node, '_sched_strategy') and node._sched_strategy is not None:
                    self.set_used_lang_feature(SchedStrategy.to_feature(node._sched_strategy))
                # Check service distribution for each class
                if hasattr(node, '_service'):
                    for jobclass, dist in node._service.items():
                        if dist is not None:
                            dist_name = type(dist).__name__
                            self.set_used_lang_feature(dist_name)
                # Check routing strategy
                if hasattr(node, '_output_strategy'):
                    for jobclass, (strategy, _, _) in node._output_strategy.items():
                        if strategy is not None:
                            self.set_used_lang_feature(RoutingStrategy.to_feature(strategy))

            elif node_type == 'Delay':
                self.set_used_lang_feature('Delay')
                self.set_used_lang_feature('SchedStrategy_INF')
                # Check service distribution
                if hasattr(node, '_service'):
                    for jobclass, dist in node._service.items():
                        if dist is not None:
                            dist_name = type(dist).__name__
                            self.set_used_lang_feature(dist_name)

            elif node_type == 'Source':
                self.set_used_lang_feature('Source')
                # Check arrival distribution
                if hasattr(node, '_arrival'):
                    for jobclass, dist in node._arrival.items():
                        if dist is not None:
                            dist_name = type(dist).__name__
                            self.set_used_lang_feature(dist_name)

            elif node_type == 'Sink':
                self.set_used_lang_feature('Sink')

            elif node_type == 'Cache':
                self.set_used_lang_feature('Cache')
                self.set_used_lang_feature('CacheClassSwitcher')
                # Check replacement strategy
                if hasattr(node, '_replacement_strategy') and node._replacement_strategy is not None:
                    self.set_used_lang_feature(ReplacementStrategy.to_feature(node._replacement_strategy))

            elif node_type == 'ClassSwitch':
                self.set_used_lang_feature('ClassSwitch')
                self.set_used_lang_feature('StatelessClassSwitcher')

            elif node_type == 'Fork':
                self.set_used_lang_feature('Fork')
                self.set_used_lang_feature('Forker')

            elif node_type == 'Join':
                self.set_used_lang_feature('Join')
                self.set_used_lang_feature('Joiner')

            elif node_type == 'Router':
                self.set_used_lang_feature('Router')
                # Check routing strategy
                if hasattr(node, '_output_strategy'):
                    for jobclass, (strategy, _, _) in node._output_strategy.items():
                        if strategy is not None:
                            self.set_used_lang_feature(RoutingStrategy.to_feature(strategy))

            elif node_type == 'Place':
                self.set_used_lang_feature('Place')
                self.set_used_lang_feature('Storage')
                self.set_used_lang_feature('Linkage')

            elif node_type == 'Transition':
                self.set_used_lang_feature('Transition')
                self.set_used_lang_feature('Enabling')
                self.set_used_lang_feature('Timing')
                self.set_used_lang_feature('Firing')

        return self._used_features

    # MATLAB-compatible aliases
    initUsedFeatures = init_used_features
    setUsedLangFeature = set_used_lang_feature
    getUsedLangFeatures = get_used_lang_features

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

    def refresh_rates(self) -> None:
        """
        Refresh only the service/arrival rates in the cached NetworkStruct.

        This is a lightweight alternative to refresh_struct() when only
        service times have changed (e.g., during iterative solver updates).
        It updates rates, scv, proc, and procid without recomputing
        routing, visits, or other structural data.

        If no cached struct exists, falls back to full refresh_struct().
        """
        if self._sn is None:
            self.refresh_struct()
            return
        self._refresh_rates()

    # CamelCase aliases
    refreshRates = refresh_rates

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
        Also clears solver results from cache nodes to prevent stale hit/miss
        probabilities from affecting subsequent solver runs.
        """
        self._has_struct = False
        self._sn = None
        self._connections = None

        # Reset cache node results to prevent stale hit/miss probabilities
        # from affecting subsequent solver runs (e.g., SSA results affecting MVA)
        from .nodes import Cache
        for node in self._nodes:
            if isinstance(node, Cache):
                node._actual_hit_prob = None
                node._actual_miss_prob = None

    def struct(self) -> 'NetworkStruct':
        """
        Get compiled NetworkStruct (alias for get_struct).

        Returns:
            NetworkStruct for solver use
        """
        return self.get_struct()

    def print_routing_matrix(self, onlyclass=None) -> None:
        """
        Print the routing matrix of the network.

        Displays routing probabilities between nodes and classes in a
        human-readable format, including class-switching routes.

        Args:
            onlyclass: Optional filter for a specific class

        References:
            MATLAB: MNetwork.printRoutingMatrix
            Java: Network.printRoutingMatrix
        """
        sn = self.get_struct()
        sn_print_routing_matrix(sn, onlyclass)

    # MATLAB-style alias
    printRoutingMatrix = print_routing_matrix

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

        # Extract class-dependent scaling
        self._refresh_class_dependence()

        # Extract joint-dependent scaling
        self._refresh_joint_dependence()

        # Extract scheduling strategies
        self._refresh_scheduling()

        # Build routing matrix
        self._refresh_routing()

        # Compute chains and visits
        self._refresh_chains()

        # Extract capacity (must be after _refresh_chains since capacity depends on chain info)
        self._refresh_capacity()

        # Extract node parameters (for transitions, caches, etc.)
        self._refresh_nodeparam()

        # Build fork-join relationship matrix
        self._refresh_fork_joins()

        # Mark struct as built before fork-join nodevisits adjustment.
        # This prevents infinite recursion: _refresh_fork_join_nodevisits()
        # calls mmt() which calls get_linked_routing_matrix() which would
        # otherwise trigger another _build_struct_from_scratch().
        self._has_struct = True

        # Adjust nodevisits for fork-join networks using MMT transformation
        self._refresh_fork_join_nodevisits()

        # Extract initial state (for SPNs)
        self._refresh_state()

        # Populate finite capacity region information
        self._refresh_regions()

    def _refresh_regions(self) -> None:
        """Populate finite capacity region information in NetworkStruct."""
        sn = self._sn
        if self._regions:
            F = len(self._regions)
            sn.nregions = F
            M = sn.nstations
            K = sn.nclasses

            sn.region = []
            sn.regionrule = np.ones((F, K), dtype=float) * DropStrategy.DROP
            sn.regionweight = np.ones((F, K), dtype=float)
            sn.regionsz = np.ones((F, K), dtype=int)

            for f, fcr in enumerate(self._regions):
                region_matrix = -1 * np.ones((M, K + 1), dtype=float)

                for node in fcr.nodes:
                    for i, station in enumerate(self._stations):
                        if station is node:
                            for r, job_class in enumerate(self._classes):
                                region_matrix[i, r] = fcr.get_class_max_jobs(job_class)
                            region_matrix[i, K] = fcr.global_max_jobs
                            break

                for r, job_class in enumerate(self._classes):
                    drop_rule = fcr.get_drop_rule(job_class)
                    sn.regionrule[f, r] = float(drop_rule)
                    sn.regionweight[f, r] = fcr.get_class_weight(job_class)
                    sn.regionsz[f, r] = fcr.get_class_size(job_class)

                sn.region.append(region_matrix)
        else:
            sn.nregions = 0
            sn.region = []
            sn.regionrule = np.array([])
            sn.regionweight = np.array([])
            sn.regionsz = np.array([])

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
                # For open classes, reference station is the Source node
                # For closed classes without explicit ref, use station 0
                is_open_class = (hasattr(jobclass, '_class_type') and
                                 jobclass._class_type == JobClassType.OPEN)
                if not is_open_class:
                    # Check if it's an OpenClass instance
                    is_open_class = jobclass.__class__.__name__ in ('OpenClass', 'OpenSignal')
                if is_open_class:
                    # Find Source station for open class
                    source_idx = 0
                    for idx, station in enumerate(self._stations):
                        if station.__class__.__name__ == 'Source':
                            source_idx = idx
                            break
                    refstat[j] = source_idx
                else:
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

        # Signal class properties
        issignal = np.zeros(nclasses, dtype=bool)
        signaltype = [None] * nclasses
        for j, jobclass in enumerate(self._classes):
            # Check if class is a Signal (or OpenSignal/ClosedSignal)
            class_name = jobclass.__class__.__name__
            if class_name in ('Signal', 'OpenSignal', 'ClosedSignal'):
                issignal[j] = True
                # Get signal type
                if hasattr(jobclass, '_signal_type') and jobclass._signal_type is not None:
                    signaltype[j] = jobclass._signal_type
                elif hasattr(jobclass, 'getSignalType'):
                    signaltype[j] = jobclass.getSignalType()
        self._sn.issignal = issignal
        self._sn.signaltype = signaltype

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
                    # Check for Immediate distribution first (mean=0)
                    # Use GlobalConstants.Immediate (large but finite value) instead of inf
                    # to match MATLAB behavior and avoid numerical issues with infinite rates
                    if hasattr(dist, 'isImmediate') and dist.isImmediate():
                        self._sn.rates[i, j] = GlobalConstants.Immediate
                    elif hasattr(dist, 'isDisabled') and dist.isDisabled():
                        # Disabled distribution: set rate=NaN to match MATLAB behavior
                        # This prevents the MVA solver from computing 1/inf=0 which differs
                        # from MATLAB's NaN handling for disabled classes
                        self._sn.rates[i, j] = np.nan
                        self._sn.scv[i, j] = np.nan
                        self._sn.procid[i, j] = ProcessType.DISABLED
                    else:
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
                    # Source, Join nodes should have rate=0 when no distribution
                    # Source: no arrival for this class
                    # Join: synchronization nodes don't do service work
                    # Cache: instant service (very high rate)
                    from .nodes import Join, Source, Cache
                    if isinstance(station, (Join, Source)):
                        self._sn.rates[i, j] = 0.0
                    elif isinstance(station, Cache):
                        # Cache has instant service for all classes
                        # Use GlobalConstants.Immediate instead of inf to match MATLAB
                        self._sn.rates[i, j] = GlobalConstants.Immediate
                        self._sn.scv[i, j] = 1.0
                    else:
                        # No service defined for this class at this station
                        # Set rate to NaN and mark as DISABLED (matching MATLAB behavior)
                        self._sn.rates[i, j] = np.nan
                        self._sn.scv[i, j] = np.nan
                        self._sn.procid[i, j] = ProcessType.DISABLED

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
            Exp, Erlang, HyperExp, Gamma, Uniform, Det, Pareto, Weibull, Lognormal, Immediate
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
        elif isinstance(dist, Immediate):
            # Check Immediate before Det since Immediate extends Det
            self._sn.procid[station_idx, class_idx] = ProcessType.IMMEDIATE
        elif isinstance(dist, Det):
            self._sn.procid[station_idx, class_idx] = ProcessType.DET
        elif isinstance(dist, Pareto):
            self._sn.procid[station_idx, class_idx] = ProcessType.PARETO
        elif isinstance(dist, Weibull):
            self._sn.procid[station_idx, class_idx] = ProcessType.WEIBULL
        elif isinstance(dist, Lognormal):
            self._sn.procid[station_idx, class_idx] = ProcessType.LOGNORMAL

        else:
            # Check for Replayer (import locally to avoid circular imports)
            from ..distributions.discrete import Replayer
            if isinstance(dist, Replayer):
                self._sn.procid[station_idx, class_idx] = ProcessType.REPLAYER
                # Store file path for JMT
                if hasattr(dist, '_file_path') and dist._file_path is not None:
                    self._sn.proc[station_idx][class_idx] = {'file_path': dist._file_path}
                return

            # Fallback: try to detect by class name
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

    def _get_process_type_id(self, dist) -> int:
        """
        Get the ProcessType ID for a distribution.

        Args:
            dist: Distribution object

        Returns:
            ProcessType ID integer
        """
        from ..distributions import (
            Exp, Erlang, HyperExp, Coxian, Det, Gamma, Uniform,
            Pareto, Weibull, Lognormal, Immediate
        )
        from ..distributions.markovian import APH, PH, MAP, MMPP2

        if dist is None:
            return ProcessType.DISABLED

        if isinstance(dist, Immediate):
            return ProcessType.IMMEDIATE
        elif isinstance(dist, Det):
            return ProcessType.DET
        elif isinstance(dist, Exp):
            return ProcessType.EXP
        elif isinstance(dist, Erlang):
            return ProcessType.ERLANG
        elif isinstance(dist, HyperExp):
            return ProcessType.HYPEREXP
        elif isinstance(dist, Coxian):
            return ProcessType.COXIAN
        elif isinstance(dist, APH):
            return ProcessType.APH
        elif isinstance(dist, PH):
            return ProcessType.PH
        elif isinstance(dist, MAP):
            return ProcessType.MAP
        elif isinstance(dist, MMPP2):
            return ProcessType.MMPP2
        elif isinstance(dist, Gamma):
            return ProcessType.GAMMA
        elif isinstance(dist, Uniform):
            return ProcessType.UNIFORM
        elif isinstance(dist, Pareto):
            return ProcessType.PARETO
        elif isinstance(dist, Weibull):
            return ProcessType.WEIBULL
        elif isinstance(dist, Lognormal):
            return ProcessType.LOGNORMAL
        else:
            # Fallback: try to detect by class name
            class_name = type(dist).__name__
            proc_type = ProcessType.fromString(class_name)
            if proc_type is not None:
                return proc_type
            return ProcessType.EXP  # Default

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
                # Handle 2D arrays: flatten or use first column
                if ld_array.ndim == 2:
                    # If 2D, use the first column (or flatten if single column/row)
                    if ld_array.shape[1] == 1:
                        ld_array = ld_array.flatten()
                    elif ld_array.shape[0] == 1:
                        ld_array = ld_array.flatten()
                    else:
                        # Multi-column: use first column as the scaling factors
                        ld_array = ld_array[:, 0]
                # Copy scaling factors, repeating last value if needed
                for j in range(total_pop):
                    if j < len(ld_array):
                        lldscaling[i, j] = ld_array[j]
                    else:
                        lldscaling[i, j] = ld_array[-1]

        self._sn.lldscaling = lldscaling

    def _refresh_class_dependence(self) -> None:
        """Extract class-dependent scaling functions from stations."""
        nstations = len(self._stations)

        # Check if any station has class dependence
        has_cd = False
        cdscaling_list = [None] * nstations

        for i, station in enumerate(self._stations):
            cd = station.get_class_dependence() if hasattr(station, 'get_class_dependence') else None
            if cd is not None:
                has_cd = True
                cdscaling_list[i] = cd

        if has_cd:
            self._sn.cdscaling = cdscaling_list
        else:
            self._sn.cdscaling = None

    def _refresh_joint_dependence(self) -> None:
        """Extract joint-dependent scaling tables and cutoffs from stations."""
        nstations = len(self._stations)
        nclasses = len(self._classes)

        has_ljd = False
        ljdscaling_list = [None] * nstations
        ljdcutoffs = np.zeros((nstations, nclasses))

        for i, station in enumerate(self._stations):
            if hasattr(station, 'get_joint_dependence'):
                scaling, cutoffs = station.get_joint_dependence()
                if scaling is not None:
                    has_ljd = True
                    ljdscaling_list[i] = np.asarray(scaling).flatten()
                    if cutoffs is not None:
                        cutoffs_arr = np.asarray(cutoffs).flatten()
                        ljdcutoffs[i, :len(cutoffs_arr)] = cutoffs_arr

        if has_ljd:
            self._sn.ljdscaling = ljdscaling_list
            self._sn.ljdcutoffs = ljdcutoffs
        else:
            self._sn.ljdscaling = None
            self._sn.ljdcutoffs = None

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
        """Extract capacity limits from stations.

        This computes capacity following MATLAB's refreshCapacity.m logic:
        - For each chain, chainCap = sum of jobs in that chain
        - classcap[ist, r] = chainCap for each class r in chain
        - Apply any explicit station/class caps
        - Final capacity = station.cap if explicit and finite,
          else min(sum(chaincap), sum(classcap))
        """
        nstations = len(self._stations)
        nclasses = len(self._classes)

        if nstations == 0 or nclasses == 0:
            self._sn.cap = np.array([])
            self._sn.classcap = np.array([]).reshape(0, 0)
            return

        # Initialize with infinity
        classcap = np.full((nstations, nclasses), np.inf)
        chaincap = np.full((nstations, nclasses), np.inf)
        capacity = np.zeros(nstations)

        # Get njobs and chains info
        njobs = self._sn.njobs.flatten() if self._sn.njobs is not None else np.zeros(nclasses)
        nchains = self._sn.nchains if hasattr(self._sn, 'nchains') and self._sn.nchains else 1
        inchain = self._sn.inchain if hasattr(self._sn, 'inchain') and self._sn.inchain else None
        rates = self._sn.rates if hasattr(self._sn, 'rates') and self._sn.rates is not None else None

        # Compute chain-based capacities
        for c in range(nchains):
            # Get classes in this chain
            if inchain is not None and c < len(inchain):
                classes_in_chain = inchain[c]
            else:
                classes_in_chain = list(range(nclasses))

            # Chain capacity = sum of jobs for all classes in chain
            # Note: Include infinite values (open classes) so that chain_cap = Inf
            # when any class has infinite jobs. This matches MATLAB's sum() behavior.
            chain_cap = np.sum([njobs[r] for r in classes_in_chain])

            for r in classes_in_chain:
                for ist in range(nstations):
                    station = self._stations[ist]

                    # Check if class is disabled at this station
                    # NaN rates means disabled, EXCEPT for Place nodes which
                    # don't have service rates (matching MATLAB refreshCapacity.m)
                    is_disabled = False
                    if rates is not None and ist < rates.shape[0] and r < rates.shape[1]:
                        if np.isnan(rates[ist, r]):
                            # Check if this is a Place node - Places are not disabled by NaN rates
                            node_idx = int(self._sn.stationToNode[ist]) if hasattr(self._sn, 'stationToNode') and self._sn.stationToNode is not None else ist
                            is_place = (node_idx < len(self._sn.nodetype) and
                                       self._sn.nodetype[node_idx] == NodeType.PLACE)
                            if not is_place:
                                is_disabled = True

                    if is_disabled:
                        classcap[ist, r] = 0
                        chaincap[ist, c] = 0
                    else:
                        chaincap[ist, c] = chain_cap
                        classcap[ist, r] = chain_cap

                        # Apply explicit class capacity if set
                        jobclass = self._classes[r] if r < len(self._classes) else None
                        if jobclass is not None:
                            class_cap = station.get_class_capacity(jobclass)
                            if np.isfinite(class_cap) and class_cap >= 0:
                                classcap[ist, r] = min(classcap[ist, r], class_cap)

                        # Apply explicit station capacity if set
                        if np.isfinite(station.capacity) and station.capacity >= 0:
                            classcap[ist, r] = min(classcap[ist, r], station.capacity)

        # Compute total capacity per station
        for ist in range(nstations):
            station = self._stations[ist]

            # If station has explicit finite cap, use it directly (Kendall K notation)
            if np.isfinite(station.capacity) and station.capacity >= 0:
                capacity[ist] = station.capacity
            else:
                # Use minimum of chain cap sum and class cap sum
                sum_chaincap = np.sum(chaincap[ist, :])
                sum_classcap = np.sum(classcap[ist, :])
                capacity[ist] = min(sum_chaincap, sum_classcap)

        # Initialize droprule matrix (default: WAITQ for all)
        # Following MATLAB's refreshCapacity.m logic:
        # - Default to WAITQ for infinite capacity
        # - Default to DROP for finite capacity (unless explicit rule set)
        droprule = np.full((nstations, nclasses), SnDropStrategy.WAITQ, dtype=np.int32)

        for ist in range(nstations):
            station = self._stations[ist]
            # Skip Source nodes - they have no drop rule
            node_idx = int(self._sn.stationToNode[ist]) if hasattr(self._sn, 'stationToNode') and self._sn.stationToNode is not None else ist
            if node_idx < len(self._sn.nodetype) and self._sn.nodetype[node_idx] == NodeType.SOURCE:
                continue

            for r in range(nclasses):
                # Check if station has explicit dropRule property
                station_drop_rule = None
                if hasattr(station, '_drop_rule'):
                    if isinstance(station._drop_rule, dict):
                        jobclass = self._classes[r] if r < len(self._classes) else None
                        if jobclass in station._drop_rule:
                            station_drop_rule = station._drop_rule[jobclass]
                    elif isinstance(station._drop_rule, list) and r < len(station._drop_rule):
                        station_drop_rule = station._drop_rule[r]

                if station_drop_rule is None:
                    # No explicit rule - default based on capacity
                    if np.isfinite(station.capacity) and station.capacity >= 0:
                        # Finite capacity: default to DROP
                        droprule[ist, r] = SnDropStrategy.DROP
                    else:
                        # Infinite capacity: keep WAITQ
                        droprule[ist, r] = SnDropStrategy.WAITQ
                else:
                    # Convert explicit rule to network_struct DropStrategy value
                    droprule[ist, r] = station_drop_rule

        self._sn.cap = capacity
        self._sn.classcap = classcap
        self._sn.droprule = droprule

    def _refresh_routing(self) -> None:
        """Build routing matrix in NetworkStruct."""
        nstations = len(self._stations)
        nclasses = len(self._classes)
        nnodes = len(self._nodes)

        if nstations == 0 or nclasses == 0:
            return

        # Transfer node._prob_routing data to _routing_matrix before building rt
        # This ensures routing probabilities set via set_prob_routing() are used
        # SKIP if _original_routes is set - this means ClassSwitch nodes were inserted
        # and _routing_matrix._routes already has the correct routes through CS nodes.
        # Adding from _prob_routing would duplicate the original routes.
        if self._routing_matrix is None or self._routing_matrix._original_routes is None:
            for node in self._nodes:
                if hasattr(node, '_prob_routing') and node._prob_routing:
                    for jobclass, routes in node._prob_routing.items():
                        for dest_node, prob in routes.items():
                            # Add to routing matrix using addRoute(class, src_node, dst_node, prob)
                            if self._routing_matrix is None:
                                self._routing_matrix = RoutingMatrix(self)
                            self._routing_matrix.addRoute(jobclass, node, dest_node, prob)

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

        # Build connection matrix early (needed for default RAND routing)
        # Use explicit links from add_link() calls as the primary source
        self._sn.connmatrix = np.zeros((nnodes, nnodes))
        for (src_idx, dst_idx) in self._links:
            if 0 <= src_idx < nnodes and 0 <= dst_idx < nnodes:
                self._sn.connmatrix[src_idx, dst_idx] = 1
        # Also add connections from routing_matrix (e.g., from set_prob_routing)
        if self._routing_matrix is not None:
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                for (node_src, node_dst), prob in routes.items():
                    if prob > 0:
                        src_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                        dst_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)
                        self._sn.connmatrix[src_idx, dst_idx] = 1

        # Build set of (node_idx, class_idx) pairs that have explicit PROB routing
        # These pairs should NOT get default RAND routing
        has_explicit_routing = set()
        # Also track (node_idx, class_idx) pairs that have INCOMING routes
        # Classes without incoming routes to a node should not have outgoing routes from it
        has_incoming_route = set()
        if self._routing_matrix is not None:
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                if isinstance(class_src, (int, np.integer)):
                    src_class_idx = class_src
                else:
                    src_class_idx = class_src._index if hasattr(class_src, '_index') else self._classes.index(class_src)
                if isinstance(class_dst, (int, np.integer)):
                    dst_class_idx = class_dst
                else:
                    dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else self._classes.index(class_dst)
                for (node_src, node_dst), prob in routes.items():
                    src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                    dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)
                    if prob > 0:
                        has_explicit_routing.add((src_node_idx, src_class_idx))
                        # Track incoming routes - destination class at destination node
                        has_incoming_route.add((dst_node_idx, dst_class_idx))

        # For closed classes, the reference station is the starting point
        # Jobs start in certain classes based on population (njobs > 0)
        # Add ONLY the reference station as "incoming" to enable routing from the starting node
        # NOTE: Previously this code added ALL Delay nodes for classes with population > 0,
        # which was wrong - jobs only start at their reference station, not all Delay nodes.
        if hasattr(self._sn, 'refstat') and self._sn.refstat is not None:
            refstat = self._sn.refstat.flatten()
            for k in range(nclasses):
                if hasattr(self._sn, 'njobs') and self._sn.njobs.size > k:
                    if self._sn.njobs.flatten()[k] > 0:
                        # Get the reference station index for this class
                        if k < len(refstat):
                            ref_station_idx = int(refstat[k])
                            # Convert station index to node index
                            if hasattr(self._sn, 'stationToNode') and self._sn.stationToNode is not None:
                                ref_node_idx = int(self._sn.stationToNode[ref_station_idx])
                            else:
                                # Fallback: use station index as node index (works for simple networks)
                                ref_node_idx = ref_station_idx
                            if ref_node_idx >= 0:
                                has_incoming_route.add((ref_node_idx, k))

        # For all classes, add all nodes with incoming connections to has_incoming_route.
        # MATLAB's getRoutingMatrix.m (lines 119-136) adds RAND routing for ALL nodes
        # that have outgoing connections, without checking if they have incoming routes.
        # We match this by adding all reachable nodes to has_incoming_route.
        # For class-switching networks, this ensures that Queue nodes downstream of
        # Delay (via RAND routing) can also route to ClassSwitch nodes.
        # NOTE: Both closed AND open classes need this - MATLAB does not discriminate.
        # Without open classes, auxiliary open classes in MMT fork-join models get
        # dead-end routing at queue nodes, producing wrong visit ratios.
        if hasattr(self._sn, 'connmatrix') and self._sn.connmatrix is not None:
            connmatrix = self._sn.connmatrix
            for ind in range(nnodes):
                # If this node has any incoming connection, add it to has_incoming_route
                # for all classes (both closed and open)
                if np.sum(connmatrix[:, ind]) > 0:
                    for k in range(nclasses):
                        has_incoming_route.add((ind, k))

        # Add default RAND routing for classes at connected nodes WITHOUT explicit routing
        # This matches MATLAB's getRoutingMatrix.m behavior at lines 119-136
        # For RAND routing (the default), MATLAB routes ALL classes from ALL
        # non-Source/non-Sink nodes to all connected destinations.
        from .nodes import Source, Sink, Cache
        from .base import RoutingStrategy

        # Find Source and Sink indices
        idx_source = None
        idx_sink = None
        for i, node in enumerate(self._nodes):
            if isinstance(node, Source):
                idx_source = i
            elif isinstance(node, Sink):
                idx_sink = i

        # Check which classes are closed (finite population)
        is_closed_class = np.isfinite(self._sn.njobs) if hasattr(self._sn, 'njobs') else np.ones(nclasses, dtype=bool)

        # Identify hit/miss classes and their source cache nodes
        # Hit/miss classes only "exist" at Cache nodes and downstream - they shouldn't
        # have routing from nodes that precede the Cache in the workflow
        hit_miss_class_indices = set()
        cache_node_indices = set()
        for i, node in enumerate(self._nodes):
            if isinstance(node, Cache):
                cache_node_indices.add(i)
                hit_classes = getattr(node, '_hit_class', {})
                miss_classes = getattr(node, '_miss_class', {})
                for hc in hit_classes.values():
                    if hc is not None:
                        hit_miss_class_indices.add(hc._index if hasattr(hc, '_index') else self._classes.index(hc))
                for mc in miss_classes.values():
                    if mc is not None:
                        hit_miss_class_indices.add(mc._index if hasattr(mc, '_index') else self._classes.index(mc))

        from .nodes import ClassSwitch as ClassSwitchNode
        for ind, node in enumerate(self._nodes):
            # Skip Source and Sink nodes for closed classes
            is_source = isinstance(node, Source)
            is_sink = isinstance(node, Sink)
            is_cache = isinstance(node, Cache)
            is_classswitch = isinstance(node, ClassSwitchNode)

            for k in range(nclasses):
                # Skip if this (node, class) has explicit routing defined (PROB routing)
                if (ind, k) in has_explicit_routing:
                    continue

                # Skip if this class has no incoming routes to this node
                # Classes that can't reach a node shouldn't have outgoing routes from it
                # Exception: Source nodes never have incoming routes but should route all open classes
                # Exception: ClassSwitch nodes need routing for all classes (Pcs matrix handles transformations)
                if (ind, k) not in has_incoming_route and not is_source and not is_classswitch:
                    continue

                # Skip hit/miss classes at non-Cache nodes (they're only created at Cache)
                # Exception: allow routing from ClassSwitch nodes that handle class switching
                # Exception: allow routing from Router nodes (they forward classes without transformation)
                if k in hit_miss_class_indices and not is_cache:
                    # Check if this is a ClassSwitch or Router node that handles this class
                    from .nodes import ClassSwitch, Router
                    if not isinstance(node, (ClassSwitch, Router)):
                        continue

                # For nodes/classes without explicit routing, add default RAND routing
                # unless node has WRROBIN with explicit weights

                # Check if this node/class uses WRROBIN with explicit weights
                jobclass = self._classes[k]
                use_wrrobin_weights = False
                wrrobin_weights = {}
                if hasattr(node, '_routing_strategies') and jobclass in node._routing_strategies:
                    from .base import RoutingStrategy as RS
                    node_strategy = node._routing_strategies[jobclass]
                    strategy_value = node_strategy.value if hasattr(node_strategy, 'value') else int(node_strategy)
                    if strategy_value == RS.WRROBIN.value:
                        # Check if we have weights for this class
                        if hasattr(node, '_routing_weights') and jobclass in node._routing_weights:
                            wrrobin_weights = node._routing_weights[jobclass]
                            if wrrobin_weights:
                                use_wrrobin_weights = True

                if is_closed_class[k]:
                    # Closed class: route from non-Source/non-Sink nodes
                    if not is_source and not is_sink:
                        # Exclude Sink from destinations for closed classes
                        connections_closed = self._sn.connmatrix[ind, :].copy()
                        if idx_sink is not None:
                            connections_closed[idx_sink] = 0
                        num_connections = np.sum(connections_closed)
                        if num_connections > 0:
                            if use_wrrobin_weights:
                                # Use WRROBIN weights to compute probabilities
                                total_weight = 0.0
                                dest_weights = {}
                                for dest_node, weight in wrrobin_weights.items():
                                    # Get destination node index
                                    dest_idx = dest_node._node_index if hasattr(dest_node, '_node_index') else self._nodes.index(dest_node)
                                    if connections_closed[dest_idx] > 0:
                                        dest_weights[dest_idx] = weight
                                        total_weight += weight
                                if total_weight > 0:
                                    for jnd, weight in dest_weights.items():
                                        self._sn.rtnodes[ind * nclasses + k, jnd * nclasses + k] = weight / total_weight
                            else:
                                # Default uniform routing
                                for jnd in range(nnodes):
                                    if connections_closed[jnd] > 0:
                                        self._sn.rtnodes[ind * nclasses + k, jnd * nclasses + k] = 1.0 / num_connections
                else:
                    # Open class: route from all connected nodes
                    num_connections = np.sum(self._sn.connmatrix[ind, :])
                    if num_connections > 0:
                        if use_wrrobin_weights:
                            # Use WRROBIN weights to compute probabilities
                            total_weight = 0.0
                            dest_weights = {}
                            for dest_node, weight in wrrobin_weights.items():
                                dest_idx = dest_node._node_index if hasattr(dest_node, '_node_index') else self._nodes.index(dest_node)
                                if self._sn.connmatrix[ind, dest_idx] > 0:
                                    dest_weights[dest_idx] = weight
                                    total_weight += weight
                            if total_weight > 0:
                                for jnd, weight in dest_weights.items():
                                    self._sn.rtnodes[ind * nclasses + k, jnd * nclasses + k] = weight / total_weight
                        else:
                            # Default uniform routing
                            for jnd in range(nnodes):
                                if self._sn.connmatrix[ind, jnd] > 0:
                                    self._sn.rtnodes[ind * nclasses + k, jnd * nclasses + k] = 1.0 / num_connections

        # Copy explicit routes from routing matrix (PROB routing)
        if self._routing_matrix is not None:
            # Copy from sparse routes to node-indexed matrix
            for (class_src, class_dst), routes in self._routing_matrix._routes.items():
                # Handle integer indices (from P[0] = ... syntax) or JobClass objects
                if isinstance(class_src, (int, np.integer)):
                    src_class_idx = class_src
                else:
                    src_class_idx = class_src._index if hasattr(class_src, '_index') else self._classes.index(class_src)
                if isinstance(class_dst, (int, np.integer)):
                    dst_class_idx = class_dst
                else:
                    dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else self._classes.index(class_dst)

                for (node_src, node_dst), prob in routes.items():
                    src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else self._nodes.index(node_src)
                    dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else self._nodes.index(node_dst)

                    i = src_node_idx * nclasses + src_class_idx
                    j = dst_node_idx * nclasses + dst_class_idx
                    self._sn.rtnodes[i, j] = prob

        # NOTE: Fork routing probabilities are left at 1.0 for each branch.
        # A Fork with arrival rate lambda sends lambda to EACH outgoing branch,
        # not lambda/fanout. The row sum > 1 is handled in sn_refresh_visits
        # which normalizes for DTMC solve and then applies fanout correction.

        # Apply ClassSwitch matrices to rtnodes
        # This matches MATLAB's getRoutingMatrix behavior for StatelessClassSwitcher
        from .nodes import ClassSwitch, Cache

        # First, apply Cache node class transformations
        # Cache nodes transform input classes to hit/miss output classes
        # Initially use 0.5/0.5 split - this will be refined during analysis
        for ind, node in enumerate(self._nodes):
            if isinstance(node, Cache):
                # Get hit/miss class mappings
                hit_classes = getattr(node, '_hit_class', {})  # Dict[input_class, output_class]
                miss_classes = getattr(node, '_miss_class', {})  # Dict[input_class, output_class]

                if hit_classes or miss_classes:
                    # Extract routing rows for this node (make a copy)
                    Pi = self._sn.rtnodes[ind * nclasses:(ind + 1) * nclasses, :].copy()

                    # Zero out routing from input classes at cache
                    # They will be replaced with routing to hit/miss classes
                    for input_class, hit_class in hit_classes.items():
                        input_idx = input_class._index if hasattr(input_class, '_index') else self._classes.index(input_class)
                        hit_idx = hit_class._index if hasattr(hit_class, '_index') else self._classes.index(hit_class)
                        miss_class = miss_classes.get(input_class)
                        miss_idx = miss_class._index if hasattr(miss_class, '_index') and miss_class else -1

                        # Get where the input class was routing to
                        for jnd in range(nnodes):
                            # Check same-class routing from input class
                            old_prob = Pi[input_idx, jnd * nclasses + input_idx]
                            if old_prob > 1e-10:
                                # Split this routing between hit and miss classes
                                # Use 0.5/0.5 split as default (will be refined in MVA)
                                self._sn.rtnodes[ind * nclasses + input_idx, jnd * nclasses + input_idx] = 0  # Remove same-class routing
                                if hit_idx >= 0:
                                    self._sn.rtnodes[ind * nclasses + input_idx, jnd * nclasses + hit_idx] = old_prob * 0.5
                                if miss_idx >= 0:
                                    self._sn.rtnodes[ind * nclasses + input_idx, jnd * nclasses + miss_idx] = old_prob * 0.5

        # Now apply ClassSwitch matrices
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

        # Route open classes back from Sink to Source
        # This creates a "closed loop" for stationary distribution computation
        # Matches MATLAB getRoutingMatrix.m lines 325-331
        has_open = any(np.isinf(self._sn.njobs))
        if has_open and idx_sink is not None and idx_source is not None:
            # Get arrival rates for open classes (from Source station)
            source_station_idx = None
            for node in self._nodes:
                if isinstance(node, Source) and hasattr(node, '_station_index'):
                    source_station_idx = node._station_index
                    break

            if source_station_idx is not None and self._sn.rates is not None:
                arv_rates = self._sn.rates[source_station_idx, :].copy()
                arv_rates[np.isnan(arv_rates)] = 0

                # Find open classes
                open_classes = np.where(np.isinf(self._sn.njobs))[0]

                # Compute chains via weakly connected components (same as MATLAB)
                # Build symmetric adjacency matrix from rtnodes
                rtnodes_sym = self._sn.rtnodes + self._sn.rtnodes.T

                # For chain detection, we use class-level connectivity
                # Two classes are in same chain if any (node,class_r) -> (node,class_s) exists
                class_adj = np.zeros((nclasses, nclasses), dtype=bool)
                for r in range(nclasses):
                    for s in range(nclasses):
                        for ind in range(nnodes):
                            for jnd in range(nnodes):
                                if rtnodes_sym[ind * nclasses + r, jnd * nclasses + s] > 0:
                                    class_adj[r, s] = True
                                    break
                            if class_adj[r, s]:
                                break

                # Find connected components using iterative approach
                visited = np.zeros(nclasses, dtype=bool)
                chains_list = []
                for start_class in range(nclasses):
                    if not visited[start_class]:
                        # BFS to find all classes in this component
                        component = []
                        queue = [start_class]
                        while queue:
                            c = queue.pop(0)
                            if not visited[c]:
                                visited[c] = True
                                component.append(c)
                                for neighbor in range(nclasses):
                                    if class_adj[c, neighbor] and not visited[neighbor]:
                                        queue.append(neighbor)
                        if component:
                            chains_list.append(component)

        # Compute stateful-indexed routing matrix (rt) from node-indexed (rtnodes)
        # using stochastic complement to absorb non-stateful nodes
        # NOTE: Compute rt BEFORE adding Sink->Source routing - solvers need the original routing.
        from ..api.mc.dtmc import dtmc_stochcomp

        # Build list of state indices corresponding to stateful nodes
        # MATLAB uses: statefulNodes = find(sn.isstateful)'
        # then builds statefulNodesClasses from those node indices
        stateful_nodes_classes = []
        if self._sn.isstateful is not None:
            stateful_node_indices = np.where(self._sn.isstateful)[0]
            for ind in stateful_node_indices:
                for k in range(nclasses):
                    stateful_nodes_classes.append(int(ind) * nclasses + k)

        if len(stateful_nodes_classes) > 0:
            stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)
            # Compute stochastic complement from rtnodes WITHOUT Sink->Source routing
            # This is used by solvers (Fluid, CTMC, etc.) for actual job routing
            rt_from_rtnodes = dtmc_stochcomp(self._sn.rtnodes, stateful_nodes_classes)
            # The stochcomp result is indexed by stateful nodes
            # rt has size (nstateful * nclasses) x (nstateful * nclasses)
            # Always use rt_from_rtnodes - it correctly handles class switching through
            # stochastic complement, absorbing inserted ClassSwitch nodes (CS_*) while
            # preserving the class transition probabilities in the routing.
            self._sn.rt = rt_from_rtnodes

        # Add Sink->Source routing to rtnodes for visit ratio calculations
        # This creates a "closed loop" making the DTMC ergodic for computing visits.
        # IMPORTANT: This is done AFTER computing rt, so solvers get the original routing.
        # We modify rtnodes in place (matching MATLAB behavior) since rt has already been computed.
        if has_open and idx_sink is not None and idx_source is not None:
            for s in open_classes:
                # Find which chain contains class s
                others_in_chain = [s]  # Default: just the class itself
                for chain in chains_list:
                    if s in chain:
                        others_in_chain = chain
                        break

                # Get arrival rates sum for classes in this chain
                chain_arv_rates = arv_rates[others_in_chain]
                total_arv = np.sum(chain_arv_rates)

                if total_arv > 0:
                    # Add routing from Sink to Source for all classes in chain
                    # rtnodes[Sink+k, Source+m] = arvRates[m] / total for all k,m in chain
                    for k in others_in_chain:
                        for m in others_in_chain:
                            prob = arv_rates[m] / total_arv if arv_rates[m] > 0 else 0
                            self._sn.rtnodes[idx_sink * nclasses + k, idx_source * nclasses + m] = prob

            # Compute rt_visits with Sink->Source for ergodic visit calculation
            if len(stateful_nodes_classes) > 0:
                rt_visits = dtmc_stochcomp(self._sn.rtnodes, stateful_nodes_classes)
                self._sn.rt_visits = rt_visits
            else:
                self._sn.rt_visits = self._sn.rt
        else:
            # For closed networks, rt and rt_visits are the same
            self._sn.rt_visits = self._sn.rt

        # Connection matrix was already built earlier in this method

        # Build routing strategy matrix (N x K) - stores routing strategy per node/class
        # Note: This is built late in the method; default RAND routing above uses RoutingStrategy.RAND
        # as the default when self._sn.routing is not yet available
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
            # If node has explicit probabilistic routing via set_prob_routing(),
            # set the routing strategy to PROB for those classes
            # This matches MATLAB behavior where setProbRouting sets RoutingStrategy.PROB
            if hasattr(node, '_prob_routing') and node._prob_routing:
                for jobclass in node._prob_routing.keys():
                    class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                    self._sn.routing[i, class_idx] = int(RoutingStrategy.PROB)

            # ClassSwitch nodes MUST use PROB routing to correctly apply class switching
            # probabilities from their switch matrix. Default RAND would produce uniform
            # class switching which is incorrect.
            from .nodes import ClassSwitch
            if isinstance(node, ClassSwitch):
                for k in range(nclasses):
                    self._sn.routing[i, k] = int(RoutingStrategy.PROB)

        # Build routing weights dictionary for WRROBIN routing
        # Dict mapping (node_idx, class_idx) -> {dest_node_idx: weight}
        self._sn.routing_weights = {}
        for i, node in enumerate(self._nodes):
            if hasattr(node, '_routing_weights') and node._routing_weights:
                for jobclass, dest_weights in node._routing_weights.items():
                    class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                    key = (i, class_idx)
                    self._sn.routing_weights[key] = {}
                    for dest_node, weight in dest_weights.items():
                        # Get destination node index
                        dest_idx = dest_node._index if hasattr(dest_node, '_index') else self._nodes.index(dest_node)
                        self._sn.routing_weights[key][dest_idx] = weight

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
        # We use rtnodes (before absorption) to capture class transitions at non-station
        # nodes like ClassSwitch, which would be lost in rt (after absorption).
        # This ensures Cache hit/miss classes are in the same chain as their input class.
        class_connected = np.zeros((nclasses, nclasses), dtype=bool)
        nnodes = len(self._nodes)
        K = nclasses

        # First, check rtnodes for class transitions (captures ClassSwitch behavior)
        if self._sn.rtnodes is not None and self._sn.rtnodes.size > 0:
            rtnodes = self._sn.rtnodes
            for r in range(K):
                for s in range(K):
                    for ind in range(nnodes):
                        for jnd in range(nnodes):
                            if rtnodes[ind * K + r, jnd * K + s] > 0:
                                class_connected[r, s] = True
                                break
                        if class_connected[r, s]:
                            break

        # Also check rt for station-level transitions (in case rtnodes isn't available)
        if self._sn.rt is not None and self._sn.rt.size > 0:
            rt = self._sn.rt
            M = nstations
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
        # Use 2D chains matrix (nchains x nclasses) for JAR/MATLAB compatibility
        # chains[chain_id, class_idx] = 1.0 if class belongs to chain, 0.0 otherwise
        chains = np.zeros((nchains, nclasses), dtype=float)
        inchain = {}

        for chain_id, root in enumerate(chain_roots):
            class_indices = chain_classes[root]
            for class_idx in class_indices:
                chains[chain_id, class_idx] = 1.0
            inchain[chain_id] = np.array(class_indices, dtype=int)

        self._sn.nchains = nchains
        self._sn.chains = chains
        self._sn.inchain = inchain

        # Force all classes in a chain to have the same reference station
        # (matches MATLAB sn_refresh_visits.m lines 48-53)
        refstat = self._sn.refstat.flatten() if self._sn.refstat is not None else np.zeros(nclasses, dtype=int)
        for chain_id in range(nchains):
            classes_in_chain = inchain[chain_id]
            if len(classes_in_chain) > 0:
                first_refstat = int(refstat[classes_in_chain[0]]) if classes_in_chain[0] < len(refstat) else 0
                for k in classes_in_chain:
                    if k < len(refstat) and int(refstat[k]) != first_refstat:
                        refstat[k] = first_refstat
        # Update sn.refstat
        if self._sn.refstat is not None:
            self._sn.refstat = refstat.reshape(self._sn.refstat.shape)
        else:
            self._sn.refstat = refstat

        # Compute visit ratios for each chain
        self._sn.visits = {}

        # Use rt_visits (stochastic complement with Sink->Source) for visits computation
        # MATLAB's sn_refresh_visits.m uses: Pchain = rt(cols,cols);
        # The rt_visits matrix is indexed by stateful nodes, with Sink→Source routing folded in
        # to make the DTMC ergodic for proper visit calculation.
        rt = self._sn.rt_visits if hasattr(self._sn, 'rt_visits') and self._sn.rt_visits is not None else self._sn.rt if self._sn.rt is not None else None

        # Update rt matrix with actual cache hit/miss probabilities if available
        # This matches MATLAB's refreshRoutingMatrix behavior (getRoutingMatrix.m lines 187-204)
        if rt is not None and rt.size > 0:
            from .base import NodeType
            for ind, node in enumerate(self._nodes):
                if hasattr(node, '_actual_hit_prob') and node._actual_hit_prob is not None:
                    # This is a cache node with actual probabilities
                    # Get the stateful index for this node
                    if hasattr(self._sn, 'nodeToStateful') and self._sn.nodeToStateful is not None:
                        isf = int(self._sn.nodeToStateful[ind])
                        if isf >= 0:
                            hit_prob = np.asarray(node._actual_hit_prob)
                            miss_prob = np.asarray(node._actual_miss_prob) if hasattr(node, '_actual_miss_prob') and node._actual_miss_prob is not None else (1.0 - hit_prob)
                            # Get hitclass/missclass arrays from nodeparam
                            hitclass = None
                            missclass = None
                            if self._sn.nodeparam is not None and ind in self._sn.nodeparam:
                                cache_param = self._sn.nodeparam[ind]
                                hitclass = np.asarray(cache_param.hitclass) if hasattr(cache_param, 'hitclass') else None
                                missclass = np.asarray(cache_param.missclass) if hasattr(cache_param, 'missclass') else None

                            if hitclass is not None and missclass is not None:
                                for r in range(nclasses):
                                    if r < len(hit_prob) and r < len(hitclass) and r < len(missclass):
                                        h = int(hitclass[r])
                                        m = int(missclass[r])
                                        if h >= 0 and m >= 0:
                                            # In the rt matrix (stochastic complement), cache class switching
                                            # is folded into routing to downstream nodes.
                                            # We need to find where rt[src_idx, :] routes to for classes h and m
                                            # and update those probabilities proportionally.
                                            src_idx = isf * nclasses + r
                                            if src_idx < rt.shape[0]:
                                                # Find all destinations for hit class h and miss class m
                                                for jsf in range(self._sn.nstateful):
                                                    hit_dst = jsf * nclasses + h
                                                    miss_dst = jsf * nclasses + m
                                                    if hit_dst < rt.shape[1] and miss_dst < rt.shape[1]:
                                                        # If there's existing routing to these destinations
                                                        if rt[src_idx, hit_dst] > 0 or rt[src_idx, miss_dst] > 0:
                                                            # Get the total probability going to this destination stateful node
                                                            old_hit = rt[src_idx, hit_dst]
                                                            old_miss = rt[src_idx, miss_dst]
                                                            old_total = old_hit + old_miss
                                                            if old_total > 0:
                                                                # Update with actual hit/miss probabilities
                                                                rt[src_idx, hit_dst] = old_total * hit_prob[r]
                                                                rt[src_idx, miss_dst] = old_total * miss_prob[r]

        if rt is not None and rt.size > 0:
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

                # Build submatrix for this chain using stateful nodes
                # MATLAB uses: Pchain = rt(cols,cols) where cols are (stateful_idx-1)*K + class_idx
                # rt is indexed by stateful nodes (M = nstateful), NOT by all nodes
                nstateful = self._sn.nstateful

                # Build indices in rt space (stateful_idx * nclasses + class_idx)
                indices = []
                for isf in range(nstateful):
                    for k in classes_in_chain:
                        indices.append(isf * nclasses + k)

                if len(indices) > 0:
                    P_chain = rt[np.ix_(indices, indices)]

                    # Use DTMC solver for both open and closed chains
                    # This matches MATLAB's sn_refresh_visits.m which uses dtmc_solve for all chains
                    from ..api.mc.dtmc import dtmc_solve, dtmc_solve_reducible
                    from .base import NodeType

                    n = len(indices)
                    row_sums = P_chain.sum(axis=1)
                    visited = row_sums > 1e-10

                    # Check for fork nodes
                    has_fork = False
                    if hasattr(self._sn, 'nodetype') and self._sn.nodetype is not None:
                        fork_type = NodeType.FORK if hasattr(NodeType, 'FORK') else getattr(NodeType, 'Fork', None)
                        if fork_type is not None:
                            # Use Python any() for enum comparison (np.any doesn't work correctly with enum lists)
                            has_fork = any(nt == fork_type for nt in self._sn.nodetype)

                    # For open chains with zero total arrival rate (e.g., auxiliary classes
                    # with Disabled arrival in MMT fork-join models), set visits to 0.
                    # In MATLAB, this produces NaN routing (0/0 in Sink→Source) which
                    # propagates through the DTMC solver, effectively giving zero visits.
                    if is_open_chain and hasattr(self._sn, 'rates') and self._sn.rates is not None:
                        from .base import NodeType as _NT
                        _source_stat = None
                        if hasattr(self._sn, 'nodetype') and self._sn.nodetype is not None:
                            for _ni in range(len(self._sn.nodetype)):
                                if self._sn.nodetype[_ni] == _NT.SOURCE:
                                    if hasattr(self._sn, 'nodeToStation') and self._sn.nodeToStation is not None:
                                        _source_stat = int(self._sn.nodeToStation[_ni])
                                    break
                        if _source_stat is not None and _source_stat >= 0:
                            _arv = self._sn.rates[_source_stat, classes_in_chain]
                            _arv = np.where(np.isnan(_arv), 0, _arv)
                            if np.sum(_arv) < 1e-10:
                                self._sn.visits[chain_id] = np.zeros((nstateful, nclasses))
                                continue

                    # Open fork networks now use the same DTMC path as closed forks.
                    # The rt_visits matrix already folds Sink→Source routing, making the
                    # DTMC ergodic. The transitive fanout correction handles fork semantics.

                    if np.sum(visited) > 0:
                        P_visited = P_chain[np.ix_(np.where(visited)[0], np.where(visited)[0])]

                        # Only normalize for Fork-containing models
                        # Fork nodes have row sums > 1 (sending to all branches with prob 1 each)
                        # which causes dtmc_solve_reducible to fail. Normalize to make stochastic.
                        # Record original row sums to apply fanout correction after DTMC solve.
                        # This matches MATLAB's sn_refresh_visits.m lines 88-98
                        row_sums_visited = np.ones(P_visited.shape[0])
                        if has_fork:
                            row_sums_visited = P_visited.sum(axis=1).flatten()
                            nonzero = (row_sums_visited > 1e-10)
                            if np.any(nonzero):
                                P_visited[nonzero] = P_visited[nonzero] / row_sums_visited[nonzero][:, np.newaxis]

                        # Try dtmc_solve first (primary), fallback to dtmc_solve_reducible
                        # This matches MATLAB's sn_refresh_visits.m line 100-106
                        # Check if chain is reducible (has transient classes)
                        from scipy.sparse.csgraph import connected_components
                        from scipy.sparse import csc_matrix
                        n_components, _ = connected_components(
                            csc_matrix(P_visited > 1e-10), directed=True, connection='strong', return_labels=True
                        )
                        is_reducible = n_components > 1

                        # Use dtmc_solve as PRIMARY, fallback to dtmc_solve_reducible
                        # This matches MATLAB's sn_refresh_visits.m lines 100-106
                        # CRITICAL: Do not change this order - it affects visit ratio computations
                        try:
                            pi_visited = dtmc_solve(P_visited)
                        except Exception:
                            pi_visited = np.full(P_visited.shape[0], np.nan)

                        # Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
                        if np.all(pi_visited == 0) or np.any(np.isnan(pi_visited)):
                            try:
                                pi_visited = dtmc_solve_reducible(P_visited)
                            except Exception:
                                pi_visited = np.ones(np.sum(visited)) / np.sum(visited)

                        pi = np.zeros(n)
                        pi[visited] = pi_visited

                        # Apply Fork fanout correction using transitive closure
                        # (matches MATLAB sn_refresh_visits.m lines 116-142)
                        # The correction must propagate through ALL transitively reachable states,
                        # not just immediate successors.
                        if has_fork:
                            # Compute transitive closure on the FULL P_chain (before visited filtering)
                            # Block outgoing edges from Join nodes to prevent the correction from
                            # leaking back through cycles (matches MATLAB sn_refresh_visits.m lines 124-133)
                            adj = (P_chain > 1e-10).astype(float)
                            join_type = NodeType.JOIN if hasattr(NodeType, 'JOIN') else getattr(NodeType, 'Join', None)
                            if join_type is not None and hasattr(self._sn, 'statefulToNode') and self._sn.statefulToNode is not None:
                                for isf in range(nstateful):
                                    nd = int(self._sn.statefulToNode[isf])
                                    if nd < len(self._sn.nodetype) and self._sn.nodetype[nd] == join_type:
                                        for k in range(n_chain_classes):
                                            adj[isf * n_chain_classes + k, :] = 0
                            reachable = adj.copy()
                            n_full = P_chain.shape[0]
                            for _ in range(int(np.ceil(np.log2(max(n_full, 2))))):
                                reachable = ((reachable + reachable @ reachable) > 0).astype(float)

                            # row_sums for the full P_chain (before normalization)
                            full_row_sums = np.array([row_sums_visited[np.searchsorted(np.where(visited)[0], row)]
                                                      if visited[row] else 1.0 for row in range(n_full)])
                            # Actually recompute from original P_chain row sums
                            # P_chain was already normalized, so we need the original row sums
                            # row_sums_visited contains the pre-normalization row sums for visited rows
                            # Map them back to full indices
                            full_row_sums = np.ones(n_full)
                            visited_idx_list = np.where(visited)[0]
                            for vi, full_idx in enumerate(visited_idx_list):
                                full_row_sums[full_idx] = row_sums_visited[vi]

                            for row in range(n_full):
                                fanout = full_row_sums[row]
                                if fanout > 1 + 1e-10:
                                    for col in range(n_full):
                                        if reachable[row, col] > 0:
                                            pi[col] *= fanout
                    else:
                        pi = np.ones(n) / n
                        row_sums_visited = np.ones(n)  # No Fork correction needed

                    # Reshape to (nstateful, nclasses) like MATLAB
                    visits_chain = np.zeros((nstateful, nclasses))
                    idx = 0
                    for isf in range(nstateful):
                        for k_idx, k in enumerate(classes_in_chain):
                            if idx < len(pi):
                                visits_chain[isf, k] = pi[idx]
                            idx += 1

                    # Normalize by SUM of visits at reference station for ALL classes in chain
                    # This matches MATLAB's sn_refresh_visits.m line 118-121:
                    # normSum = sum(visits{c}(sn.stationToStateful(refstat(inchain{c}(1))),inchain{c}))
                    ref_class = classes_in_chain[0]
                    ref_stat = int(self._sn.refstat[ref_class])  # station index
                    ref_sf = int(self._sn.stationToStateful[ref_stat])  # stateful index
                    norm_sum = np.sum(visits_chain[ref_sf, classes_in_chain])
                    if norm_sum > 1e-10:
                        visits_chain = visits_chain / norm_sum
                    elif is_reducible and has_fork and is_open_chain:
                        # For open fork networks with absorbing states (Sink),
                        # compute visits from the routing structure.
                        # Each forked branch has visits = 1/fanout.
                        visits_chain = self._compute_fork_visits(
                            nstateful, nclasses, classes_in_chain, ref_sf
                        )
                    elif is_reducible:
                        # Reference station has 0 steady-state visits (transient state).
                        # For reducible chains with absorbing states, MATLAB falls back
                        # to uniform visits (all stations have visit ratio = 1).
                        # This matches MATLAB's behavior which warns about reducibility
                        # but still computes results using uniform visit ratios.
                        import warnings
                        warnings.warn(
                            "Reducible network topology detected, results may be unreliable.",
                            UserWarning
                        )
                        # Set uniform visits for all stations in this chain
                        for isf in range(nstateful):
                            for k in classes_in_chain:
                                visits_chain[isf, k] = 1.0

                    self._sn.visits[chain_id] = visits_chain
                else:
                    self._sn.visits[chain_id] = np.ones((nstateful, nclasses))
        else:
            # No routing - default visits
            nstateful = self._sn.nstateful if hasattr(self._sn, 'nstateful') else nstations
            for chain_id in range(nchains):
                self._sn.visits[chain_id] = np.ones((nstateful, nclasses))

        # Compute node visits (for non-station nodes like Router, Sink)
        # This is needed for computing arrival rates at non-station nodes
        self._sn.nodevisits = {}
        nnodes = len(self._nodes)

        if self._sn.rtnodes is not None and self._sn.rtnodes.size > 0:
            from ..api.mc.dtmc import dtmc_solve, dtmc_solve_reducible

            for chain_id in range(nchains):
                if chain_id not in self._sn.inchain:
                    continue

                classes_in_chain = self._sn.inchain[chain_id].flatten().astype(int)
                n_chain_classes = len(classes_in_chain)

                # Build indices for nodes in this chain
                nodes_cols = []
                for ind in range(nnodes):
                    for ik, k in enumerate(classes_in_chain):
                        nodes_cols.append(ind * nclasses + k)

                nodes_cols = np.array(nodes_cols, dtype=int)

                # Extract routing submatrix for this chain
                if len(nodes_cols) > 0:
                    nodes_Pchain = self._sn.rtnodes[np.ix_(nodes_cols, nodes_cols)].copy()

                    # Handle NaN values in routing matrix
                    for row in range(nodes_Pchain.shape[0]):
                        nan_cols = np.isnan(nodes_Pchain[row, :])
                        if np.any(nan_cols):
                            non_nan_sum = np.sum(nodes_Pchain[row, ~nan_cols])
                            remaining_prob = max(0, 1 - non_nan_sum)
                            n_nan = np.sum(nan_cols)
                            if n_nan > 0 and remaining_prob > 0:
                                nodes_Pchain[row, nan_cols] = remaining_prob / n_nan
                            else:
                                nodes_Pchain[row, nan_cols] = 0

                    # Find visited nodes
                    nodes_visited = np.sum(nodes_Pchain, axis=1) > 1e-10

                    if np.sum(nodes_visited) > 0:
                        nodes_Pchain_visited = nodes_Pchain[np.ix_(np.where(nodes_visited)[0],
                                                                   np.where(nodes_visited)[0])]

                        # Normalize rows, recording original row sums for fork correction
                        # (matches MATLAB sn_refresh_visits.m lines 198-209)
                        nodes_row_sums_visited = nodes_Pchain_visited.sum(axis=1).flatten()
                        nonzero = (nodes_row_sums_visited > 1e-10)
                        if np.any(nonzero):
                            nodes_Pchain_visited[nonzero] = nodes_Pchain_visited[nonzero] / nodes_row_sums_visited[nonzero][:, np.newaxis]

                        # Solve DTMC for node visits
                        try:
                            nodes_alpha_visited = dtmc_solve(nodes_Pchain_visited)
                        except Exception:
                            nodes_alpha_visited = np.zeros(nodes_Pchain_visited.shape[0])

                        if np.all(nodes_alpha_visited == 0) or np.any(np.isnan(nodes_alpha_visited)):
                            try:
                                nodes_alpha_visited = dtmc_solve_reducible(nodes_Pchain_visited)
                            except Exception:
                                nodes_alpha_visited = np.ones(np.sum(nodes_visited)) / np.sum(nodes_visited)

                        nodes_alpha = np.zeros(len(nodes_cols))
                        nodes_alpha[nodes_visited] = nodes_alpha_visited

                        # Apply Fork fanout correction for node visits
                        # Block propagation through Join nodes (matches MATLAB sn_refresh_visits.m lines 220-252)
                        if has_fork:
                            n_nodes_full = nodes_Pchain.shape[0]
                            nodes_adj = (nodes_Pchain > 1e-10).astype(float)
                            join_type = NodeType.JOIN if hasattr(NodeType, 'JOIN') else getattr(NodeType, 'Join', None)
                            if join_type is not None and hasattr(self._sn, 'nodetype') and self._sn.nodetype is not None:
                                for nd in range(nnodes):
                                    if nd < len(self._sn.nodetype) and self._sn.nodetype[nd] == join_type:
                                        for k in range(n_chain_classes):
                                            nodes_adj[nd * n_chain_classes + k, :] = 0
                            nodes_reachable = nodes_adj.copy()
                            for _ in range(int(np.ceil(np.log2(max(n_nodes_full, 2))))):
                                nodes_reachable = ((nodes_reachable + nodes_reachable @ nodes_reachable) > 0).astype(float)

                            # Map row sums back to full node indices
                            nodes_full_row_sums = np.ones(n_nodes_full)
                            nodes_visited_idx_list = np.where(nodes_visited)[0]
                            for vi, full_idx in enumerate(nodes_visited_idx_list):
                                nodes_full_row_sums[full_idx] = nodes_row_sums_visited[vi]

                            for row in range(n_nodes_full):
                                fanout = nodes_full_row_sums[row]
                                if fanout > 1 + 1e-10:
                                    for col in range(n_nodes_full):
                                        if nodes_reachable[row, col] > 0:
                                            nodes_alpha[col] *= fanout
                    else:
                        nodes_alpha = np.zeros(len(nodes_cols))

                    # Reshape to (nnodes, nclasses)
                    nodevisits_chain = np.zeros((nnodes, nclasses))
                    idx = 0
                    for ind in range(nnodes):
                        for k in classes_in_chain:
                            if idx < len(nodes_alpha):
                                nodevisits_chain[ind, k] = nodes_alpha[idx]
                            idx += 1

                    # Normalize by reference station visits
                    ref_class = classes_in_chain[0]
                    ref_stat = int(self._sn.refstat[ref_class])
                    ref_node = self._sn.stationToNode[ref_stat] if ref_stat < len(self._sn.stationToNode) else 0
                    node_norm_sum = np.sum(nodevisits_chain[ref_node, classes_in_chain])
                    if node_norm_sum > 1e-10:
                        nodevisits_chain = nodevisits_chain / node_norm_sum

                    nodevisits_chain[nodevisits_chain < 0] = 0
                    nodevisits_chain[np.isnan(nodevisits_chain)] = 0

                    self._sn.nodevisits[chain_id] = nodevisits_chain

        # Set reference class per chain
        # In MATLAB, refclass=0 means "no reference class" (1-indexed arrays)
        # In Python (0-indexed), we use -1 to mean "no reference class"
        # When refclass >= 0, sn_get_residt_from_respt uses that class's visits
        # at the reference station as the divisor for WN computation.
        # MATLAB sets refclass by finding classes with isrefclass=True for each chain.
        refclass = -np.ones(nchains, dtype=int)

        # Find reference classes from is_reference_class attribute
        refclasses = np.zeros(nclasses, dtype=bool)
        for k, cls in enumerate(self._classes):
            if hasattr(cls, 'is_reference_class') and cls.is_reference_class:
                refclasses[k] = True
            elif hasattr(cls, '_is_reference_class') and cls._is_reference_class:
                refclasses[k] = True

        # For each chain, find the reference class (intersection of inchain and refclasses)
        for c in range(nchains):
            if self._sn.inchain is not None and c in self._sn.inchain:
                inchain_classes = self._sn.inchain[c]
                refclass_in_chain = np.intersect1d(inchain_classes, np.where(refclasses)[0])
                if len(refclass_in_chain) > 0:
                    # Use the first reference class found in this chain
                    refclass[c] = refclass_in_chain[0]

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

        # Handle Queue nodes with polling and switchover
        from .nodes import Queue
        from ..api.sn import SchedStrategy

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Queue):
                # Check for polling type or switchover settings
                polling_type = node.get_polling_type() if hasattr(node, 'get_polling_type') else None
                switchover = getattr(node, '_switchover', None)

                if polling_type is not None or switchover:
                    # Initialize per-class parameters
                    if node_idx not in nodeparam:
                        nodeparam[node_idx] = {}

                    for r, jobclass in enumerate(self._classes):
                        if r not in nodeparam[node_idx]:
                            nodeparam[node_idx][r] = {}

                        # Store polling type and parameters
                        if polling_type is not None:
                            nodeparam[node_idx][r]['pollingType'] = polling_type
                            nodeparam[node_idx][r]['pollingPar'] = [getattr(node, '_polling_k', 1)]

                        # Store switchover times
                        if switchover:
                            # Check for class-to-class switchover (tuple key) or single-class (class key)
                            switchover_times = {}
                            switchover_proc_ids = {}

                            for key, dist in switchover.items():
                                if isinstance(key, tuple):
                                    # (from_class, to_class) -> distribution
                                    from_class, to_class = key
                                    if from_class == jobclass:
                                        to_idx = self._classes.index(to_class) if to_class in self._classes else -1
                                        if to_idx >= 0:
                                            switchover_times[to_idx] = dist
                                            switchover_proc_ids[to_idx] = self._get_process_type_id(dist)
                                elif key == jobclass:
                                    # Single class switchover (for polling)
                                    switchover_times[0] = dist
                                    switchover_proc_ids[0] = self._get_process_type_id(dist)

                            if switchover_times:
                                nodeparam[node_idx][r]['switchoverTime'] = switchover_times
                                nodeparam[node_idx][r]['switchoverProcId'] = switchover_proc_ids

        # Handle Cache nodes
        from .nodes import Cache

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Cache):
                @dataclass
                class CacheParam:
                    """Container for cache parameters."""
                    nitems: int = 0
                    cap: int = 0
                    itemcap: np.ndarray = field(default_factory=lambda: np.array([]))
                    hitclass: np.ndarray = field(default_factory=lambda: np.array([]))
                    missclass: np.ndarray = field(default_factory=lambda: np.array([]))
                    actualhitprob: Optional[np.ndarray] = None
                    actualmissprob: Optional[np.ndarray] = None
                    replacestrat: Any = None  # Renamed from 'replacement' to match MATLAB
                    accost: Optional[np.ndarray] = None
                    pread: List[Any] = field(default_factory=list)  # PMF values per class

                param = CacheParam()
                param.nitems = node._num_items if hasattr(node, '_num_items') else 0
                # _item_level_cap is a list/array - store as itemcap and take first as cap
                if hasattr(node, '_item_level_cap') and node._item_level_cap is not None:
                    item_cap = node._item_level_cap
                    if isinstance(item_cap, (list, np.ndarray)) and len(item_cap) > 0:
                        param.itemcap = np.asarray(item_cap)
                        param.cap = int(item_cap[0]) if not hasattr(item_cap[0], 'item') else int(item_cap[0])
                    else:
                        param.itemcap = np.array([item_cap])
                        param.cap = int(item_cap)
                else:
                    param.itemcap = np.array([0])
                    param.cap = 0
                param.replacestrat = node._replacement_strategy if hasattr(node, '_replacement_strategy') else None

                # Populate pread - PMF values from read distribution for each class
                param.pread = [None] * nclasses
                if hasattr(node, '_read_process') and node._read_process:
                    for jobclass, dist in node._read_process.items():
                        if jobclass in self._classes and dist is not None:
                            class_idx = self._classes.index(jobclass)
                            # Evaluate PMF for items 1 to nitems
                            if hasattr(dist, 'evalPMF'):
                                pmf_values = np.array([dist.evalPMF(i) for i in range(1, param.nitems + 1)])
                                param.pread[class_idx] = pmf_values
                            elif hasattr(dist, 'pmf'):
                                # Try scipy-style pmf method
                                pmf_values = np.array([dist.pmf(i) for i in range(1, param.nitems + 1)])
                                param.pread[class_idx] = pmf_values

                # Initialize hitclass and missclass arrays
                hitclass = np.zeros(nclasses, dtype=int) - 1  # -1 means no hit class
                missclass = np.zeros(nclasses, dtype=int) - 1  # -1 means no miss class

                # Get hit/miss class mappings from the cache node
                if hasattr(node, '_hit_class') and node._hit_class:
                    for in_class, out_class in node._hit_class.items():
                        in_idx = self._classes.index(in_class) if in_class in self._classes else -1
                        out_idx = self._classes.index(out_class) if out_class in self._classes else -1
                        if in_idx >= 0 and out_idx >= 0:
                            hitclass[in_idx] = out_idx
                if hasattr(node, '_miss_class') and node._miss_class:
                    for in_class, out_class in node._miss_class.items():
                        in_idx = self._classes.index(in_class) if in_class in self._classes else -1
                        out_idx = self._classes.index(out_class) if out_class in self._classes else -1
                        if in_idx >= 0 and out_idx >= 0:
                            missclass[in_idx] = out_idx

                param.hitclass = hitclass
                param.missclass = missclass

                # Access cost (if available)
                if hasattr(node, '_accost') and node._accost is not None:
                    param.accost = np.array(node._accost)

                nodeparam[node_idx] = param

        # Handle Logger nodes
        from .nodes import Logger

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Logger):
                @dataclass
                class LoggerParam:
                    """Container for logger parameters."""
                    fileName: str = 'log.csv'
                    filePath: str = '/tmp/'
                    startTime: bool = False
                    loggerName: bool = False
                    timestamp: bool = True
                    jobID: bool = True
                    jobClass: bool = True
                    timeSameClass: bool = False
                    timeAnyClass: bool = False

                param = LoggerParam()
                param.fileName = node.file_name if hasattr(node, 'file_name') else 'log.csv'
                param.filePath = node.file_path if hasattr(node, 'file_path') else '/tmp/'
                param.startTime = node._want_start_time if hasattr(node, '_want_start_time') else False
                param.loggerName = node._want_logger_name if hasattr(node, '_want_logger_name') else False
                param.timestamp = node._want_timestamp if hasattr(node, '_want_timestamp') else True
                param.jobID = node._want_job_id if hasattr(node, '_want_job_id') else True
                param.jobClass = node._want_job_class if hasattr(node, '_want_job_class') else True
                param.timeSameClass = node._want_time_same_class if hasattr(node, '_want_time_same_class') else False
                param.timeAnyClass = node._want_time_any_class if hasattr(node, '_want_time_any_class') else False

                nodeparam[node_idx] = param

        # Handle Fork nodes with tasks per link (fanOut)
        from .nodes import Fork
        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Fork):
                tasks_per_link = node.get_tasks_per_link()
                if tasks_per_link is not None:
                    # Convert to numpy array for consistent handling
                    tasks_per_link = np.atleast_1d(tasks_per_link)
                    # Use the first value (or max if multiple values exist)
                    fanout_val = int(tasks_per_link.flat[0]) if tasks_per_link.size > 0 else 1
                else:
                    fanout_val = 1
                nodeparam[node_idx] = {'fanOut': fanout_val}

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

    def _refresh_fork_join_nodevisits(self) -> None:
        """
        Adjust nodevisits for fork-join networks using MMT transformation.

        In MATLAB, this is done inline in refreshStruct.m (lines 423-462).
        For networks with fork-join nodes, the nodevisits are recomputed
        using the mixed-model transformation (MMT) from ModelAdapter.

        The algorithm:
        1. Call ModelAdapter.mmt() to get transformed model without forks
        2. For each new chain in the transformed model (auxiliary chains):
           - Find the original fork and chain
           - Get auxiliary class visits from transformed model
           - Add scaled visits to original class visits
        """
        sn = self._sn

        # Check if there are any fork-join relationships
        if sn.fj is None or not np.any(sn.fj):
            return

        # Import ModelAdapter for MMT transformation
        try:
            from ..io.model_adapter import ModelAdapter
        except ImportError:
            return

        # Perform MMT transformation
        try:
            mmt_result = ModelAdapter.mmt(self)
        except Exception:
            # MMT failed - skip adjustment
            return

        nonfjmodel = mmt_result.nonfjmodel
        fjclassmap = mmt_result.fjclassmap
        forkmap = mmt_result.fjforkmap
        fanout = mmt_result.fanout

        # If any fanout is 1, warn about partial support (matches MATLAB line 429-433)
        if np.any(fanout == 1):
            import warnings
            warnings.warn(
                "The specified fork-join topology has partial support, "
                "only SolverJMT simulation results may be reliable.",
                UserWarning
            )

        # Get struct from transformed model
        nonfjmodel.refresh_struct()
        fsn = nonfjmodel._sn

        if fsn is None:
            return

        # Process each new chain (auxiliary chains created by MMT)
        # In Python, fjclassmap/forkmap are 0-indexed arrays where index i
        # corresponds to auxiliary class at index (sn.nclasses + i) in fsn
        for aux_idx in range(len(fjclassmap)):
            # Find which chain this auxiliary class belongs to in fsn
            aux_class = sn.nclasses + aux_idx
            new_chain = None
            for c in range(fsn.nchains):
                if c in fsn.inchain and aux_class in fsn.inchain[c]:
                    new_chain = c
                    break

            if new_chain is None:
                continue

            # Get original fork and chain
            orig_fork = int(forkmap[aux_idx])
            orig_class = int(fjclassmap[aux_idx])

            # Find original chain containing orig_class
            orig_chain = None
            for c in range(sn.nchains):
                if c in sn.inchain and orig_class in sn.inchain[c]:
                    orig_chain = c
                    break

            if orig_chain is None:
                continue

            if new_chain not in fsn.nodevisits:
                continue

            # Get auxiliary visits and zero out Source/Sink/Fork (matches MATLAB line 440)
            Vaux = fsn.nodevisits[new_chain].copy()
            from .base import NodeType
            for i, nt in enumerate(fsn.nodetype):
                if nt == NodeType.SOURCE or nt == NodeType.SINK or nt == NodeType.FORK:
                    Vaux[i, :] = 0

            # Get only columns for classes in the auxiliary chain
            aux_classes = list(fsn.inchain[new_chain])
            Vaux = Vaux[:, aux_classes]

            # If node counts differ, map Vaux by node name (matches MATLAB lines 442-456)
            if fsn.nnodes != sn.nnodes:
                VauxMapped = np.zeros((sn.nnodes, Vaux.shape[1]))
                for fsn_row in range(fsn.nnodes):
                    fsn_node_name = fsn.nodenames[fsn_row] if fsn.nodenames else None
                    if fsn_node_name is None:
                        continue
                    for sn_row in range(sn.nnodes):
                        sn_node_name = sn.nodenames[sn_row] if sn.nodenames else None
                        if sn_node_name == fsn_node_name:
                            VauxMapped[sn_row, :] = Vaux[fsn_row, :]
                            break
                Vaux = VauxMapped

            if orig_chain not in sn.nodevisits:
                continue

            # Get original visits
            X = sn.nodevisits[orig_chain].copy()

            # Get fanOut from nodeparam (matches MATLAB line 460)
            fanOut_val = 1
            if sn.nodeparam is not None and orig_fork < len(sn.nodeparam):
                if sn.nodeparam[orig_fork] is not None:
                    np_fanOut = getattr(sn.nodeparam[orig_fork], 'fanOut', None)
                    if np_fanOut is not None and np_fanOut > 0:
                        fanOut_val = np_fanOut

            # Add scaled auxiliary visits to original (matches MATLAB lines 458-460)
            orig_classes = list(sn.inchain[orig_chain])
            for jaux, j in enumerate(orig_classes):
                if jaux < Vaux.shape[1]:
                    # Get fanout for this auxiliary class
                    aux_fanout = fanout[aux_idx] if aux_idx < len(fanout) else 1
                    if aux_fanout == 0 or np.isnan(aux_fanout):
                        aux_fanout = 1
                    scaled_Vaux = Vaux[:, jaux] * aux_fanout
                    self._sn.nodevisits[orig_chain][:, j] = fanOut_val * (X[:, j] + scaled_Vaux)

    def _compute_fork_visits(
        self, nstateful: int, nclasses: int, classes_in_chain: np.ndarray, ref_sf: int,
        rt: np.ndarray = None
    ) -> np.ndarray:
        """
        Compute visits for open fork networks with absorbing states.

        In a fork network, a Fork with arrival rate lambda sends lambda to
        EACH outgoing branch, not lambda/fanout. Therefore, the visits at
        each branch station should be 1.0 (same as the reference station).

        Note: Fork nodes are typically not stateful, so the fork behavior is
        represented in the routing matrix as row sums > 1 (e.g., Source routes
        to multiple queues each with probability 1).

        Args:
            nstateful: Number of stateful nodes
            nclasses: Number of classes
            classes_in_chain: Array of class indices in this chain
            ref_sf: Stateful index of reference station (typically Source)
            rt: Routing matrix (optional, defaults to self._sn.rt)

        Returns:
            visits_chain: Array of shape (nstateful, nclasses) with visits
        """
        visits_chain = np.zeros((nstateful, nclasses))

        # The routing matrix rt is indexed by (stateful * nclasses + class)
        # Fork behavior is represented as stations with row sums > 1
        if rt is None:
            rt = self._sn.rt
        if rt is None:
            # Fallback: set uniform visits
            for isf in range(nstateful):
                for k in classes_in_chain:
                    visits_chain[isf, k] = 1.0
            return visits_chain

        # Find stations that act as fork sources (row sum > 1 for any class in chain)
        # and their successors (fork branches)
        fork_sources = {}  # fork_source_sf -> {fanout, successors}

        for src_sf in range(nstateful):
            successors = set()
            for k in classes_in_chain:
                src_idx = src_sf * nclasses + k
                if src_idx < rt.shape[0]:
                    # Find all destinations with non-zero probability
                    for dest_sf in range(nstateful):
                        for dest_k in classes_in_chain:
                            dest_idx = dest_sf * nclasses + dest_k
                            if dest_idx < rt.shape[1] and rt[src_idx, dest_idx] > 1e-10:
                                if dest_sf != src_sf:
                                    successors.add(dest_sf)

            # Check if this is a fork source (more than 1 successor)
            # For fork networks, the row sum indicates simultaneous routing to all branches
            row_sum = 0
            for k in classes_in_chain:
                src_idx = src_sf * nclasses + k
                if src_idx < rt.shape[0]:
                    row_sum = max(row_sum, np.sum(rt[src_idx, :]))

            if len(successors) > 1 and row_sum > 1 + 1e-10:
                # This is a fork source - routes to multiple destinations simultaneously
                fork_sources[src_sf] = {
                    'fanout': len(successors),
                    'successors': list(successors)
                }

        # If no fork sources found, set uniform visits
        if not fork_sources:
            for isf in range(nstateful):
                for k in classes_in_chain:
                    visits_chain[isf, k] = 1.0
            return visits_chain

        # Set visits for each stateful node
        # - Reference station (Source): visits = 1
        # - Stations that are fork branches: visits = 1 (Fork sends full rate to each branch)
        # - Other stations: visits = 1

        for isf in range(nstateful):
            for k in classes_in_chain:
                if isf == ref_sf:
                    # Reference station (Source) - check if class actually arrives here
                    # For class switching in open networks, only the originating class
                    # has arrivals at the Source; switched-to classes have visits = 0
                    ref_stat = int(self._sn.statefulToStation[ref_sf]) if hasattr(self._sn, 'statefulToStation') else ref_sf
                    if (hasattr(self._sn, 'rates') and self._sn.rates is not None and
                        ref_stat < self._sn.rates.shape[0] and k < self._sn.rates.shape[1]):
                        arrival_rate = self._sn.rates[ref_stat, k]
                        if np.isnan(arrival_rate) or arrival_rate <= 0:
                            # No arrivals for this class at Source
                            visits_chain[isf, k] = 0.0
                        else:
                            visits_chain[isf, k] = 1.0
                    else:
                        visits_chain[isf, k] = 1.0
                else:
                    # Check if this station is a successor of any fork source
                    is_fork_branch = False
                    total_fanout = 1
                    for fork_sf, fork_info in fork_sources.items():
                        if isf in fork_info['successors']:
                            is_fork_branch = True
                            total_fanout = fork_info['fanout']
                            break

                    if is_fork_branch:
                        # Fork branch: visits = 1 (Fork sends full arrival rate to each branch)
                        visits_chain[isf, k] = 1.0
                    else:
                        # Other stations: visits = 1
                        visits_chain[isf, k] = 1.0

        return visits_chain

    def _refresh_state(self) -> None:
        """
        Extract initial state from stateful nodes.

        This populates sn.state with initial token counts for Places in SPNs.
        For closed classes, jobs start at their reference stations.
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

        # Track which classes have explicit state set
        explicit_state_set = np.zeros(nclasses, dtype=bool)

        # Track which stateful nodes had their state explicitly set
        explicit_node_set = set()

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, StatefulNode):
                # Get stateful index
                stateful_idx = None
                if nodeToStateful is not None and node_idx < len(nodeToStateful):
                    stateful_idx = int(nodeToStateful[node_idx])

                if stateful_idx is not None and stateful_idx >= 0 and stateful_idx < len(state):
                    # Check if this node had setState() called explicitly
                    explicitly_set = getattr(node, '_state_explicitly_set', False)
                    # Get node state
                    node_state = node.get_state() if hasattr(node, 'get_state') else node._state if hasattr(node, '_state') else None
                    if node_state is not None:
                        node_state = np.asarray(node_state).flatten()
                        if explicitly_set:
                            # Use the full state as-is (including zeros)
                            for k in range(min(nclasses, len(node_state))):
                                state[stateful_idx][k] = node_state[k]
                                explicit_state_set[k] = True
                            explicit_node_set.add(stateful_idx)
                        else:
                            # Only copy non-zero values (legacy behavior)
                            for k in range(min(nclasses, len(node_state))):
                                if node_state[k] > 0:
                                    state[stateful_idx][k] = node_state[k]
                                    explicit_state_set[k] = True

        # For closed classes without explicit state, place jobs at reference station
        stationToStateful = self._sn.stationToStateful if hasattr(self._sn, 'stationToStateful') else None
        if stationToStateful is not None:
            stationToStateful = np.asarray(stationToStateful).flatten()

        refstat = self._sn.refstat if hasattr(self._sn, 'refstat') else None
        if refstat is not None:
            refstat = np.asarray(refstat).flatten()

        njobs = self._sn.njobs if hasattr(self._sn, 'njobs') else None
        if njobs is not None:
            njobs = np.asarray(njobs).flatten()

        if refstat is not None and njobs is not None and stationToStateful is not None:
            for r in range(nclasses):
                if not explicit_state_set[r] and np.isfinite(njobs[r]) and njobs[r] > 0:
                    # Closed class without explicit state - place at reference station
                    ref_station = int(refstat[r])
                    if ref_station < len(stationToStateful):
                        stateful_idx = int(stationToStateful[ref_station])
                        if stateful_idx >= 0 and stateful_idx < len(state):
                            state[stateful_idx][r] = njobs[r]

        # Override with stored marginals if available (from initDefault/initFromMarginal).
        # This ensures sn.state[isf] always contains correct per-class marginal counts,
        # even for FCFS queues whose raw state encodes buffer class IDs + phases.
        # Matches MATLAB where sn.state{isf} marginals are computed via State.toMarginal
        # inside getStruct() -> getState().
        state_marginal = getattr(self, '_state_marginal', None)
        if state_marginal is not None and stationToStateful is not None:
            nstations = len(self._stations) if hasattr(self, '_stations') else 0
            marginal_2d = state_marginal.reshape(nstations, nclasses) if len(state_marginal) == nstations * nclasses else None
            if marginal_2d is not None:
                for ist in range(nstations):
                    isf = int(stationToStateful[ist])
                    if isf >= 0 and isf < len(state):
                        state[isf] = marginal_2d[ist, :].copy()

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

    def _state_from_marginal_and_started(self, station, n_vec, s_vec) -> np.ndarray:
        """
        Generate state vector for a station from marginal queue lengths and started jobs.

        The state format depends on the scheduling strategy:
        - FCFS/HOL/LCFS: Ordered buffer with class IDs for each job position
        - PS/INF/DPS/GPS: Job counts per class
        - SIRO: Unordered buffer counts + service state

        Args:
            station: The station node
            n_vec: Vector of queue lengths per class at this station
            s_vec: Vector of jobs in service per class at this station

        Returns:
            State vector for the station
        """
        nclasses = len(self._classes)
        n_vec = np.asarray(n_vec).flatten()
        s_vec = np.asarray(s_vec).flatten()

        # Get scheduling strategy value (use .value for comparison to handle
        # different SchedStrategy enum imports from constants vs base modules)
        sched_raw = getattr(station, '_sched_strategy', SchedStrategy.INF)
        sched = sched_raw.value if hasattr(sched_raw, 'value') else sched_raw

        # Define scheduling strategy groups by their integer values
        # FCFS=0, LCFS=1, LCFSPR=2, LCFSPI=3, HOL=9
        fcfs_group = {SchedStrategy.FCFS.value, SchedStrategy.HOL.value,
                      SchedStrategy.LCFS.value, SchedStrategy.LCFSPR.value,
                      SchedStrategy.LCFSPI.value}
        # INF=7, PS=4, DPS=5, GPS=6, LPS=17, PSPRIO=21, DPSPRIO=19, GPSPRIO=20
        inf_group = {SchedStrategy.INF.value, SchedStrategy.PS.value,
                     SchedStrategy.DPS.value, SchedStrategy.GPS.value,
                     SchedStrategy.LPS.value, SchedStrategy.PSPRIO.value,
                     SchedStrategy.DPSPRIO.value, SchedStrategy.GPSPRIO.value}
        # SIRO=12, SEPT=10, LEPT=11, POLLING=15
        siro_group = {SchedStrategy.SIRO.value, SchedStrategy.SEPT.value,
                      SchedStrategy.LEPT.value, SchedStrategy.POLLING.value}

        # FCFS, HOL, LCFS: ordered buffer representation
        # Each job is represented by its class ID (1-indexed for compatibility with MATLAB)
        if sched in fcfs_group:
            total_jobs = int(np.sum(n_vec))
            if total_jobs == 0:
                # Empty state: just zeros for phases
                return np.zeros(max(nclasses, 1))

            # Build buffer: class IDs (1-indexed) for each job in the queue
            # Jobs in service are represented at the end, buffer jobs first
            buffer = []
            for r in range(nclasses):
                njobs_class = int(n_vec[r])
                # Add class ID (1-indexed) for each job of this class
                for _ in range(njobs_class):
                    buffer.append(r + 1)  # 1-indexed class ID

            return np.array(buffer, dtype=float)

        # PS, INF (Delay), DPS, GPS, LPS, PSPRIO, DPSPRIO, GPSPRIO: count per class
        elif sched in inf_group:
            return n_vec.copy()

        # SIRO, SEPT, LEPT, etc.: unordered buffer + service state
        elif sched in siro_group:
            # [buffer counts per class, service phase per class]
            return n_vec.copy()

        # Default: just return counts
        else:
            return n_vec.copy()

    def init_default(self) -> None:
        """
        Initialize network with default state.

        For closed classes, all jobs start at their reference station.
        For open classes, sources start with potential arrivals.

        Delegates to initFromMarginal to generate proper per-node state spaces
        and state priors, which are needed for CTMC transient analysis.
        """
        nclasses = len(self._classes)
        nstations = len(self._stations)

        # Compute initial marginal: all closed class jobs at their reference station
        n0 = np.zeros((nstations, nclasses))

        for r, jobclass in enumerate(self._classes):
            # Check if this is a closed class with a reference station
            if hasattr(jobclass, '_refstat') and jobclass._refstat is not None:
                refstat = jobclass._refstat
                # Find the station index
                try:
                    ist = self._stations.index(refstat)
                    njobs = getattr(jobclass, '_njobs', 0)
                    if np.isfinite(njobs):
                        n0[ist, r] = njobs
                except ValueError:
                    pass  # Station not found

        # Delegate to initFromMarginal which generates per-node state spaces
        # and state priors (needed for CTMC/FLD transient analysis)
        self.init_from_marginal(n0)

    def has_init_state(self) -> bool:
        """Check if network has initialized state."""
        return any(node.state is not None for node in self._nodes if node.is_stateful())

    def get_state(self) -> list:
        """
        Get current state of all stateful nodes.

        Returns:
            List of state arrays, one per stateful node, in node order.
            Each element is a list containing the node's state array(s).
        """
        if not self.has_init_state():
            self.init_default()
        state = []
        for node in self._nodes:
            if node.is_stateful():
                s = node.get_state()
                state.append([s] if s is not None else [None])
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

    def get_state_marginal(self) -> np.ndarray:
        """
        Get the stored marginal state used for CTMC initial state matching.

        Returns:
            Flat array of marginal job counts [n[0,0], n[0,1], ..., n[M-1,K-1]]
            where n[i,k] is number of jobs of class k at station i.
            Returns None if no marginal has been set.
        """
        return getattr(self, '_state_marginal', None)

    # PascalCase alias
    getStateMarginal = get_state_marginal

    def init_from_marginal_and_started(self, n, s) -> None:
        """
        Initialize network state from marginal queue lengths and started jobs.

        This method generates the full state space for each station based on
        the marginal job counts, sets the state prior (first state with
        probability 1 by default), and sets the current state.

        Args:
            n: Marginal queue lengths matrix (nstations x nclasses).
               n[i][r] is the number of jobs of class r at station i.
            s: Started jobs matrix (nstations x nclasses).
               s[i][r] is the number of jobs of class r in service at station i.
        """
        from ..api.state.marginal import fromMarginal, fromMarginalAndStarted

        n = np.atleast_2d(n)
        s = np.atleast_2d(s)

        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Validate dimensions
        if n.shape[0] != nstations:
            raise ValueError(f"n must have {nstations} rows (one per station)")
        if s.shape[0] != nstations:
            raise ValueError(f"s must have {nstations} rows (one per station)")

        # Get network struct for state space generation
        sn = self.getStruct()

        # Store the requested marginal for later use (e.g., CTMC solver initial state)
        # Flatten n into a vector [n[0,0], n[0,1], ..., n[1,0], n[1,1], ...] = M*K elements
        # Pad n to nclasses columns if needed
        n_padded = np.zeros((nstations, nclasses))
        for i in range(min(n.shape[0], nstations)):
            for k in range(min(n.shape[1], nclasses)):
                n_padded[i, k] = n[i, k]
        self._state_marginal = n_padded.flatten()

        # Set state for each station using scheduling-aware state generation
        for i, station in enumerate(self._stations):
            if hasattr(station, 'set_state'):
                n_row = n[i, :] if n.shape[1] >= nclasses else np.zeros(nclasses)
                s_row = s[i, :] if s.shape[1] >= nclasses else np.zeros(nclasses)
                for k in range(min(nclasses, n.shape[1])):
                    n_row[k] = n[i, k]
                for k in range(min(nclasses, s.shape[1])):
                    s_row[k] = s[i, k]

                # Get node index for this station
                node_idx = self._nodes.index(station) if station in self._nodes else i

                # Generate state space for this station
                # Use fromMarginalAndStarted if started jobs specified, otherwise fromMarginal
                if np.any(s_row > 0):
                    state_space = fromMarginalAndStarted(sn, node_idx, n_row, s_row)
                else:
                    state_space = fromMarginal(sn, node_idx, n_row)

                if state_space is not None and len(state_space) > 0:
                    # Set state space on the station
                    if hasattr(station, 'set_state_space'):
                        station.set_state_space(state_space)

                    # Set state prior: first state has probability 1
                    if hasattr(station, 'setStatePrior'):
                        if len(state_space) == 1:
                            station.setStatePrior(np.array([1.0]))
                        else:
                            # Multiple states: first state gets probability 1
                            prior = np.zeros(len(state_space))
                            prior[0] = 1.0
                            station.setStatePrior(prior)

                    # Set current state to first state in space
                    station.set_state(state_space[0])
                else:
                    # Fallback to simple state generation
                    state_vec = self._state_from_marginal_and_started(station, n_row, s_row)
                    station.set_state(state_vec)

        self._has_state = True
        self._reset_struct()

    # PascalCase alias
    initFromMarginalAndStarted = init_from_marginal_and_started

    def init_from_marginal(self, n) -> None:
        """
        Initialize network state from marginal queue lengths only.

        This is equivalent to calling init_from_marginal_and_started with
        a zero matrix for started jobs.

        Args:
            n: Marginal queue lengths. Can be:
               - 1D list/array with one value per station (single class)
               - 2D list/array (nstations x nclasses)
        """
        n = np.atleast_2d(n)

        # Handle 1D input (single class - vector of length nstations)
        if n.shape[0] == 1 and len(self._stations) > 1:
            # If shape is (1, nstations), transpose to (nstations, 1)
            n = n.T

        # Create zero started matrix
        s = np.zeros_like(n)

        self.init_from_marginal_and_started(n, s)

    # PascalCase alias
    initFromMarginal = init_from_marginal

    # PascalCase alias
    initDefault = init_default
    getState = get_state
    setState = set_state

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
                # After deepcopy, _fork is a copy of the old node, not the old node itself.
                # Match by node index to find the correct node in the new network.
                fork_idx = new_node._fork._node_index if hasattr(new_node._fork, '_node_index') else -1
                if 0 <= fork_idx < len(new_network._nodes):
                    new_node._fork = new_network._nodes[fork_idx]
                elif new_node._fork in node_map:
                    new_node._fork = node_map[new_node._fork]

            # Update service/arrival process class references
            # Note: After deepcopy, the keys are new class objects, so we need to
            # match by index to find the corresponding new class in the network
            if hasattr(new_node, '_service_process') and new_node._service_process:
                new_service = {}
                for deepcopied_cls, dist in new_node._service_process.items():
                    # Find the matching new class by index
                    cls_idx = deepcopied_cls._index if hasattr(deepcopied_cls, '_index') else None
                    if cls_idx is not None and cls_idx < len(new_network._classes):
                        new_service[new_network._classes[cls_idx]] = dist
                    else:
                        new_service[deepcopied_cls] = dist
                new_node._service_process = new_service

            if hasattr(new_node, '_arrival_process') and new_node._arrival_process:
                new_arrival = {}
                for deepcopied_cls, dist in new_node._arrival_process.items():
                    # Find the matching new class by index
                    cls_idx = deepcopied_cls._index if hasattr(deepcopied_cls, '_index') else None
                    if cls_idx is not None and cls_idx < len(new_network._classes):
                        new_arrival[new_network._classes[cls_idx]] = dist
                    else:
                        new_arrival[deepcopied_cls] = dist
                new_node._arrival_process = new_arrival

            # Update routing strategies class references
            if hasattr(new_node, '_routing_strategies') and new_node._routing_strategies:
                new_strategies = {}
                for deepcopied_cls, strategy in new_node._routing_strategies.items():
                    # Find the matching new class by index
                    cls_idx = deepcopied_cls._index if hasattr(deepcopied_cls, '_index') else None
                    if cls_idx is not None and cls_idx < len(new_network._classes):
                        new_strategies[new_network._classes[cls_idx]] = strategy
                    else:
                        new_strategies[deepcopied_cls] = strategy
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

        # Copy explicit links (these are node index tuples, so just copy the set)
        new_network._links = set(self._links)

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
    # TRANSIENT METRIC HANDLES
    # =====================================================================

    def getTranHandles(self):
        """Get transient metric handles.

        Returns:
            Tuple (Qt, Ut, Tt) of nested lists of Metric objects for
            transient queue-length, utilization, and throughput.
        """
        from .nodes import Source, Sink, Fork, Join
        from ..constants import Metric, MetricType

        M = self.get_number_of_stations()
        K = self.get_number_of_classes()

        Qt = [[None for _ in range(K)] for _ in range(M)]
        Ut = [[None for _ in range(K)] for _ in range(M)]
        Tt = [[None for _ in range(K)] for _ in range(M)]

        for ist in range(M):
            station = self._stations[ist]
            for r in range(K):
                job_class = self._classes[r]
                Qt[ist][r] = Metric(MetricType.TranQLen, job_class, station)
                Ut[ist][r] = Metric(MetricType.TranUtil, job_class, station)
                Tt[ist][r] = Metric(MetricType.TranTput, job_class, station)
                if isinstance(station, (Source, Sink)):
                    Qt[ist][r].disabled = True
                    Ut[ist][r].disabled = True
                if isinstance(station, (Fork, Join)):
                    Ut[ist][r].disabled = True

        return Qt, Ut, Tt

    get_tran_handles = getTranHandles

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
            probability of routing from node i to node j. The matrix is sized
            for all nodes in the model, with indices matching model.get_nodes() order.
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

        # Try to get the model from the first node to determine full node list
        # This ensures the matrix is sized correctly for the network
        model = None
        for node in nodes:
            if hasattr(node, '_model') and node._model is not None:
                model = node._model
                break

        if model is not None:
            # Get all nodes from model in their canonical order
            all_nodes = model.get_nodes()
            n = len(all_nodes)
            # Build mapping from node object to its index in model
            node_to_model_idx = {node: idx for idx, node in enumerate(all_nodes)}
        else:
            # Fallback: use only the nodes in the path (legacy behavior)
            unique_nodes = []
            node_to_model_idx = {}
            for node in nodes:
                if node not in node_to_model_idx:
                    node_to_model_idx[node] = len(unique_nodes)
                    unique_nodes.append(node)
            n = len(unique_nodes)

        # Create 2D list of routing probabilities
        P = [[0.0] * n for _ in range(n)]

        # Forward routing: node[i] -> node[i+1]
        for i in range(len(nodes) - 1):
            from_node = nodes[i]
            to_node = nodes[i + 1]
            if from_node in node_to_model_idx and to_node in node_to_model_idx:
                from_idx = node_to_model_idx[from_node]
                to_idx = node_to_model_idx[to_node]
                P[from_idx][to_idx] = 1.0

        # Close the loop for non-sink networks (closed networks)
        # Only if the last node is not a Sink
        from .nodes import Sink
        if len(nodes) > 0 and not isinstance(nodes[-1], Sink):
            last_node = nodes[-1]
            first_node = nodes[0]
            if last_node in node_to_model_idx and first_node in node_to_model_idx:
                last_idx = node_to_model_idx[last_node]
                first_idx = node_to_model_idx[first_node]
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
                node = Delay(model, f'Station{i+1}')
            else:
                node = Queue(model, f'Station{i+1}', strategy)
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
