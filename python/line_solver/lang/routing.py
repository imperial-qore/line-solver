"""
Routing matrix for LINE queueing network models (pure Python).

This module provides the RoutingMatrix class for defining how jobs
move between nodes in a stochastic network.
"""

from typing import Any, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
import numpy as np

from .base import JobClass, Node

if TYPE_CHECKING:
    from .network import Network


class RoutingMatrix:
    """
    Matrix representing routing probabilities between network nodes.

    The routing matrix defines how jobs move between nodes in the network,
    specifying the probability that a job leaving one node will arrive at
    another node.

    This is a pure Python implementation that stores routing data in numpy
    arrays and can be converted to Java when needed.

    Args:
        network: The parent network this routing matrix belongs to.
    """

    def __init__(self, network: 'Network'):
        """
        Initialize a routing matrix for a network.

        Args:
            network: The parent network.
        """
        self.network = network
        self._nodes_cache = None

        # Store routing as sparse structure: (class_src, class_dst) -> {(node_src, node_dst): prob}
        self._routes: Dict[Tuple[JobClass, JobClass], Dict[Tuple[Node, Node], float]] = {}

        # Store original routes before ClassSwitch node insertion
        # Used by toMatrix() to build the rt matrix with correct class-switching probabilities
        self._original_routes: Optional[Dict[Tuple[JobClass, JobClass], Dict[Tuple[Node, Node], float]]] = None

        # Also maintain dense matrix for efficiency
        self._matrix: Optional[np.ndarray] = None
        self._java_obj = None

    def _get_node_by_index(self, index: Union[int, Node]) -> Node:
        """
        Convert node index to node object.

        Args:
            index: Node index or node object.

        Returns:
            Node object corresponding to the index.
        """
        if self._nodes_cache is None:
            self._nodes_cache = self.network.get_nodes()

        if isinstance(index, (int, np.integer)):
            if 0 <= index < len(self._nodes_cache):
                return self._nodes_cache[index]
            else:
                raise IndexError(f"Node index {index} out of range (0-{len(self._nodes_cache)-1})")
        return index

    def set(self, *args) -> 'RoutingMatrix':
        """
        Set routing probabilities in the matrix.

        Supports multiple argument patterns:
        - 5 args: set(class_source, class_dest, node_source, node_dest, probability)
        - 3 args: set(class_source, class_dest, routing_matrix)
        - 2 args: set(jobclass, routing_matrix)

        Returns:
            Self for method chaining.
        """
        if len(args) == 5:
            class_source, class_dest, node_source, node_dest, prob = args
            node_source = self._get_node_by_index(node_source)
            node_dest = self._get_node_by_index(node_dest)

            key = (class_source, class_dest)
            if key not in self._routes:
                self._routes[key] = {}
            self._routes[key][(node_source, node_dest)] = prob

        elif len(args) == 3:
            class_source, class_dest, rt = args
            # Convert list to numpy array if needed
            if isinstance(rt, list):
                rt = np.asarray(rt)
            if isinstance(rt, np.ndarray):
                nodes = self.network.get_nodes()
                for i in range(rt.shape[0]):
                    for j in range(rt.shape[1]):
                        if rt[i, j] != 0:
                            self.set(class_source, class_dest, nodes[i], nodes[j], rt[i, j])
            elif isinstance(rt, RoutingMatrix):
                # Copy routes from another matrix
                for (cs, cd), routes in rt._routes.items():
                    if cs == class_source and cd == class_dest:
                        for (ns, nd), prob in routes.items():
                            self.set(class_source, class_dest, ns, nd, prob)

        elif len(args) == 2:
            jobclass, rt = args
            if isinstance(rt, list):
                rt = np.asarray(rt)
            if isinstance(rt, np.ndarray):
                self.set(jobclass, jobclass, rt)
            elif isinstance(rt, RoutingMatrix):
                self.set(jobclass, jobclass, rt)
        else:
            raise ValueError(f"Unsupported number of arguments: {len(args)}. Expected 2, 3, or 5.")

        self._matrix = None  # Invalidate cached matrix
        self._java_obj = None
        return self

    def __getitem__(self, index: int) -> 'RoutingMatrixRowView':
        """
        Enable [i][j] access to routing probabilities by node index.

        Args:
            index: Row index (source node index).

        Returns:
            RoutingMatrixRowView for second-level indexing.
        """
        return RoutingMatrixRowView(self, index)

    def _resolve_class(self, key):
        """Resolve an integer index or JobClass to a JobClass object."""
        if isinstance(key, (int, np.integer)):
            classes = self.network.get_classes()
            if 0 <= key < len(classes):
                return classes[key]
            raise IndexError(f"Class index {key} out of range (0-{len(classes)-1})")
        return key

    def __setitem__(self, key, value):
        """
        Set routing using indexing notation.

        Args:
            key: Single jobclass (or int index) or tuple of (source_class, dest_class).
            value: 2D numpy array of routing probabilities.
        """
        if isinstance(key, tuple) and len(key) == 1:
            source_class = dest_class = self._resolve_class(key[0])
        elif isinstance(key, tuple) and len(key) == 2:
            source_class = self._resolve_class(key[0])
            dest_class = self._resolve_class(key[1])
        elif not isinstance(key, tuple):
            source_class = dest_class = self._resolve_class(key)
        else:
            raise ValueError("Key must be a single jobclass or tuple of (source_class, dest_class)")

        # Handle RoutingMatrix objects directly
        if isinstance(value, RoutingMatrix):
            self.set(source_class, dest_class, value)
            return

        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) != 2:
            raise ValueError("Value must be a 2D array representing routing probabilities")

        nodes = self.network.get_nodes()
        if value.shape[0] != len(nodes) or value.shape[1] != len(nodes):
            raise ValueError(f"Routing matrix must be {len(nodes)}x{len(nodes)} to match network topology")

        for i in range(value.shape[0]):
            for j in range(value.shape[1]):
                if value[i, j] != 0:
                    self.set(source_class, dest_class, nodes[i], nodes[j], float(value[i, j]))

    def addRoute(self, jobclass: JobClass, *args):
        """
        Add a routing path through multiple nodes for a job class.

        Args:
            jobclass: Job class to configure routing for.
            *args: Nodes followed by optional probability.
                   If last arg is a number, it's used as probability.
        """
        if len(args) < 2:
            raise ValueError("addRoute requires at least 2 nodes (source and destination)")

        # Check if last arg is probability
        if len(args) >= 2 and isinstance(args[-1], (int, float, np.integer, np.floating)) and not hasattr(args[-1], '_node_index'):
            nodes = args[:-1]
            probability = args[-1]
        else:
            nodes = args
            probability = 1.0

        if len(nodes) < 2:
            raise ValueError("addRoute requires at least 2 nodes (source and destination)")

        # Create chain of routes
        for i in range(len(nodes) - 1):
            self.set(jobclass, jobclass, nodes[i], nodes[i + 1], probability)

    def addClassSwitch(self, source_class: JobClass, dest_class: JobClass,
                       source_node: Node, dest_node: Node, probability: float):
        """
        Add a class switching route between nodes.

        Args:
            source_class: The job class before the switch.
            dest_class: The job class after the switch.
            source_node: The node where the job departs from.
            dest_node: The node where the job arrives.
            probability: The probability of taking this route.
        """
        self.set(source_class, dest_class, source_node, dest_node, probability)

    def setRoutingMatrix(self, jobclass: Union[JobClass, List[JobClass]],
                         nodes: List[Node], pmatrix: np.ndarray):
        """
        Set routing probabilities using a matrix for specific job class(es).

        Args:
            jobclass: Job class or list of job classes.
            nodes: List of nodes in the routing matrix.
            pmatrix: 2D or 3D matrix of routing probabilities.
        """
        if isinstance(jobclass, JobClass):
            # 2D matrix for single class
            for i in range(len(nodes)):
                for j in range(len(nodes)):
                    if pmatrix[i][j] != 0:
                        self.set(jobclass, jobclass, nodes[i], nodes[j], pmatrix[i][j])
        else:
            # 3D matrix for multiple classes
            for k, jc in enumerate(jobclass):
                for i in range(len(nodes)):
                    for j in range(len(nodes)):
                        if pmatrix[k][i][j] != 0:
                            self.set(jc, jc, nodes[i], nodes[j], pmatrix[k][i][j])

    def toMatrix(self) -> np.ndarray:
        """
        Convert the routing to a dense matrix.

        Non-station nodes (like Router, ClassSwitch) are absorbed by computing
        transitive routing through them.

        Returns:
            ndarray: (M*K) x (M*K) routing matrix indexed by station then class.
        """
        if self._matrix is not None:
            return self._matrix

        nstations = self.network.get_number_of_stations()
        nclasses = self.network.get_number_of_classes()
        nnodes = self.network.get_number_of_nodes()

        if nstations == 0 or nclasses == 0:
            return np.zeros((0, 0))

        nodes = self.network.get_nodes()
        classes = self.network.get_classes()

        # Build node to station mapping and identify station/non-station nodes
        node_to_station = {}
        station_nodes = []
        non_station_nodes = []
        for i, node in enumerate(nodes):
            if hasattr(node, '_station_index') and node._station_index is not None and node._station_index >= 0:
                node_to_station[node] = node._station_index
                station_nodes.append(i)
            else:
                non_station_nodes.append(i)

        # Build station-class indexed matrix
        self._matrix = np.zeros((nstations * nclasses, nstations * nclasses))

        # Use original routes if available (before ClassSwitch node insertion)
        # This preserves the correct class-switching probabilities
        routes_to_use = self._original_routes if self._original_routes is not None else self._routes

        # For each class pair, build node-level routing and absorb non-stations
        for (class_src, class_dst), routes in routes_to_use.items():
            # Handle integer indices (from P[0] = ... syntax) or JobClass objects
            if isinstance(class_src, (int, np.integer)):
                src_class_idx = class_src
            else:
                src_class_idx = class_src._index if hasattr(class_src, '_index') else classes.index(class_src)
            if isinstance(class_dst, (int, np.integer)):
                dst_class_idx = class_dst
            else:
                dst_class_idx = class_dst._index if hasattr(class_dst, '_index') else classes.index(class_dst)

            # Build node-level routing matrix for this class pair
            P_nodes = np.zeros((nnodes, nnodes))
            for (node_src, node_dst), prob in routes.items():
                # Handle integer indices (from P[0] = ... syntax) or Node objects
                if isinstance(node_src, (int, np.integer)):
                    src_node_idx = node_src
                else:
                    src_node_idx = node_src._node_index if hasattr(node_src, '_node_index') else nodes.index(node_src)
                if isinstance(node_dst, (int, np.integer)):
                    dst_node_idx = node_dst
                else:
                    dst_node_idx = node_dst._node_index if hasattr(node_dst, '_node_index') else nodes.index(node_dst)
                P_nodes[src_node_idx, dst_node_idx] = prob

            # Absorb non-station nodes
            # For each non-station node n, update: P[i,j] += P[i,n] * P[n,j] / (1 - P[n,n])
            # Do this iteratively until all non-station routing is absorbed
            for _ in range(len(non_station_nodes) + 1):  # Multiple passes for chains of non-stations
                for n in non_station_nodes:
                    # Get routing through this non-station node
                    p_in = P_nodes[:, n].copy()  # Probability of reaching node n
                    p_out = P_nodes[n, :].copy()  # Probability of leaving node n
                    p_self = P_nodes[n, n]  # Self-loop probability

                    if np.sum(p_in) < 1e-10 or np.sum(p_out) < 1e-10:
                        continue

                    # Compute effective routing through n (accounting for self-loops)
                    scale = 1.0 / (1.0 - p_self) if p_self < 1.0 - 1e-10 else 1.0

                    # Add transitive routing: for all i,j, P[i,j] += P[i,n] * P[n,j] * scale
                    for i in range(nnodes):
                        if p_in[i] > 1e-10 and i != n:
                            for j in range(nnodes):
                                if p_out[j] > 1e-10 and j != n:
                                    P_nodes[i, j] += p_in[i] * p_out[j] * scale

                    # Zero out the non-station node's routing
                    P_nodes[:, n] = 0
                    P_nodes[n, :] = 0

            # Extract station-to-station routing
            for src_node in station_nodes:
                for dst_node in station_nodes:
                    prob = P_nodes[src_node, dst_node]
                    if prob > 1e-10:
                        src_station = node_to_station[nodes[src_node]]
                        dst_station = node_to_station[nodes[dst_node]]
                        i = src_station * nclasses + src_class_idx
                        j = dst_station * nclasses + dst_class_idx
                        self._matrix[i, j] = prob

        return self._matrix

    # Aliases
    set_routing_matrix = setRoutingMatrix
    add_route = addRoute
    add_class_switch = addClassSwitch


class RoutingMatrixRowView:
    """
    Helper class to enable [i][j] access on RoutingMatrix.

    This allows accessing routing probabilities with matrix[i][j] notation
    where i is the source node index and j is the destination node index.
    """

    def __init__(self, routing_matrix: RoutingMatrix, row_index: int):
        self._routing_matrix = routing_matrix
        self._row_index = row_index

    def __getitem__(self, col_index: int) -> float:
        """
        Get routing probability from row node to column node.

        Args:
            col_index: Column index (destination node index).

        Returns:
            Routing probability (summed across all class pairs).
        """
        rm = self._routing_matrix
        nodes = rm.network.get_nodes()

        if self._row_index < 0 or self._row_index >= len(nodes):
            raise IndexError(f"Row index {self._row_index} out of range")
        if col_index < 0 or col_index >= len(nodes):
            raise IndexError(f"Column index {col_index} out of range")

        src_node = nodes[self._row_index]
        dst_node = nodes[col_index]

        # Sum probabilities across all class pairs
        total_prob = 0.0
        for routes in rm._routes.values():
            if (src_node, dst_node) in routes:
                total_prob += routes[(src_node, dst_node)]

        return total_prob

    def __setitem__(self, col_index: int, value) -> None:
        """
        Set routing for a class pair using P[class_src_idx][class_dst_idx] = matrix syntax.

        This enables MATLAB-style routing matrix assignment where P{i,j} = matrix
        is translated to Python as P[i][j] = matrix.

        Args:
            col_index: Destination class index (0-based).
            value: 2D list or numpy array of routing probabilities between nodes.
        """
        rm = self._routing_matrix
        classes = rm.network.get_classes()

        if self._row_index < 0 or self._row_index >= len(classes):
            raise IndexError(f"Source class index {self._row_index} out of range (0-{len(classes)-1})")
        if col_index < 0 or col_index >= len(classes):
            raise IndexError(f"Destination class index {col_index} out of range (0-{len(classes)-1})")

        src_class = classes[self._row_index]
        dst_class = classes[col_index]

        # Convert to numpy array if needed
        if isinstance(value, list):
            value = np.array(value)

        if isinstance(value, np.ndarray):
            rm.set(src_class, dst_class, value)
        elif isinstance(value, RoutingMatrix):
            rm.set(src_class, dst_class, value)
        else:
            raise TypeError(f"Unsupported value type: {type(value)}. Expected list, numpy array, or RoutingMatrix.")
