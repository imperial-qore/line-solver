"""
Network class for LINE native Python implementation.

This module provides the Network class for constructing and managing queueing networks.
Ported from MATLAB implementation in matlab/src/lang/@MNetwork/
"""

from typing import Optional, Dict, List, Tuple, Union
import numpy as np
import scipy.sparse as sp

from .base import (
    Element, ElementType, Network as NetworkBase, NetworkElement, Node, StatefulNode, Station,
    JobClass, JobClassType, NodeType, SchedStrategy, RoutingStrategy, DropStrategy
)
from .region import Region
from .routing import RoutingMatrix
from ..api.sn.network_struct import NetworkStruct, DropStrategy as SnDropStrategy
from ..api.sn.utils import sn_print_routing_matrix
from ..api.io import console as _console
from ..constants import ProcessType, GlobalConstants


def lld_scaling_as_1d(ld):
    """Reduce a user-supplied lldScaling to the 1-D alpha(n) the struct holds.

    setLoadDependence takes a vector indexed by the station population, but a
    caller may hand in a column, a row, or a (n, R) table with one column per
    class. Everything that reads the scaling -- sn.lldscaling, the JSON wire
    form -- has to reduce it the SAME way, or a model solves natively and
    fails once it crosses a codebase boundary.
    """
    ld_array = np.asarray(ld)
    if ld_array.ndim == 2:
        # single row/column flattens; multi-column takes the first column
        if ld_array.shape[1] == 1 or ld_array.shape[0] == 1:
            ld_array = ld_array.flatten()
        else:
            ld_array = ld_array[:, 0]
    return ld_array


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
        self._gd_scaling = None  # global (Whittle) rate scaling phi(n); see set_global_dependence
        self._gd_scaling_cutoff = 10  # open-class wire truncation of _gd_scaling; solving ignores it
        self._gd_scaling_peak = None
        self._has_struct = False  # Whether struct is valid
        self._rates_dirty = False  # Whether service/arrival rates need refresh
        self._rewards = {}  # Dict[reward_name, reward_function]
        self._source_idx = -1  # Cached source node index
        self._sink_idx = -1  # Cached sink node index
        self._log_path = None  # Path for logger output files
        self.allow_replace = False  # When True, add_node replaces nodes with same name
        self._do_checks = True  # When False, skip link-time model checks (setChecks)
        self.fork_stateful = False  # Fork nodes are stateful (FJ-augmented copies only, see ModelAdapter.fjtag)
        self.is_fj_augmented = False  # True on fork-join tag-augmented copies (skips the MMT visit correction)

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

    def remove_class(self, jobclass) -> None:
        """
        Remove a job class from the network, in place.

        Every node drops the per-class configuration referencing the class (see
        Node.remove_job_class), the class is dropped from the model, the
        remaining classes are re-indexed and the struct is invalidated so that
        the next solve recompiles it. Mirrors MATLAB @MNetwork/removeClass.m and
        Network.removeClass in the JAR. Use ModelAdapter.remove_class (or
        Network.without_class) for the non-mutating variant.

        Args:
            jobclass: JobClass to remove, or its 0-based index

        Raises:
            RuntimeError: if the class is the only class of the network, or the
                model contains a Cache node
        """
        if isinstance(jobclass, (int, np.integer)):
            idx = int(jobclass)
            if idx < 0 or idx >= len(self._classes):
                raise ValueError(
                    f"Invalid class index {idx}, must be in [0, {len(self._classes) - 1}]")
            target = self._classes[idx]
        else:
            target = next((c for c in self._classes if c is jobclass), None)
            if target is None:
                # a class object of a copied model resolves by name, matching
                # @MNetwork/removeClass.m
                name = getattr(jobclass, 'name', None)
                target = next((c for c in self._classes if c.name == name), None)

        if len(self._classes) <= 1:
            if target is not None:
                raise RuntimeError(
                    "The network has a single class, it cannot be removed from the model.")
            return
        if target is None:
            return

        for node in self._nodes:
            node.remove_job_class(target)
        if getattr(self, '_routing_matrix', None) is not None:
            self._routing_matrix.remove_job_class(target)
        self._classes.remove(target)
        for i, jclass in enumerate(self._classes):
            jclass._set_index(i)
        self.reset()

    removeClass = remove_class

    def without_class(self, jobclass) -> 'Network':
        """
        Return a copy of this network with a job class removed.

        The original model is left untouched, so it stays solvable; use this for
        ablation studies or a per-class decomposition. Mirrors
        ModelAdapter.removeClass in MATLAB and the JAR.

        Args:
            jobclass: JobClass to remove, or its 0-based index

        Returns:
            Network: a copy without the specified class
        """
        newmodel = self.copy()
        newmodel.remove_class(jobclass)
        return newmodel

    withoutClass = without_class

    def aggregate_chains(self, suffix: str = ''):
        """
        Return a copy of this model with one aggregate class per chain.

        All classes belonging to the same chain are merged into a single class,
        which eliminates class switching. The aggregation is exact for
        product-form models and approximate otherwise, since one aggregate
        service process replaces the per-class ones weighted by alpha.

        Args:
            suffix: suffix appended to the aggregate class names

        Returns:
            tuple: (chain_model, alpha, deagg_info), where alpha is the
                (nstations, nclasses) matrix of aggregation factors and
                deagg_info carries the data needed by
                sn_deaggregate_chain_results to map chain-level metrics back to
                class-level metrics
        """
        from ..api.io.model_adapter import ModelAdapter
        return ModelAdapter.aggregate_chains(self, suffix)

    aggregateChains = aggregate_chains

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

        # tracks connectivity only; routing probabilities computed later in _refresh_routing (mirrors MATLAB addLink).

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

        # Flatten list/tuple arguments (MATLAB-style addRegion({n1,n2}) compat)
        flat_nodes = []
        for nd in node_list:
            if isinstance(nd, (list, tuple)):
                flat_nodes.extend(nd)
            else:
                flat_nodes.append(nd)
        node_list = flat_nodes

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

    def set_checks(self, val: bool) -> None:
        """Enable or disable link-time model checks (e.g. the OI permutation-
        invariance check). Mirrors MATLAB model.setChecks."""
        self._do_checks = bool(val)

    setChecks = set_checks

    def _iter_routing_blocks(self, routing_matrix):
        """
        Normalize the accepted ``link`` argument shapes to ``(r, s, arr)``
        blocks, where ``r``/``s`` are class indices and ``arr`` is an
        (nnodes x nnodes) float array of the routing probabilities from class
        ``r`` to class ``s``.

        Yielding blocks rather than re-sniffing the format at each call site
        keeps the validation below in step with link's own conversion branches.
        """
        nodes = self.get_nodes()
        classes = self.get_classes()
        nnodes = len(nodes)
        nclasses = len(classes)

        def _as2d(m):
            if m is None:
                return None
            try:
                arr = np.asarray(m, dtype=float)
            except (ValueError, TypeError):
                return None
            return arr if arr.ndim == 2 and arr.size else None

        if isinstance(routing_matrix, RoutingMatrix):
            cidx = {c: i for i, c in enumerate(classes)}
            nidx = {n: i for i, n in enumerate(nodes)}
            blocks = {}
            for (csrc, cdst), routes in routing_matrix._routes.items():
                r, s = cidx.get(csrc), cidx.get(cdst)
                if r is None or s is None:
                    continue
                arr = blocks.setdefault((r, s), np.zeros((nnodes, nnodes)))
                for (nsrc, ndst), prob in routes.items():
                    i, j = nidx.get(nsrc), nidx.get(ndst)
                    if i is not None and j is not None:
                        arr[i, j] = prob
            for (r, s), arr in blocks.items():
                yield r, s, arr

        elif isinstance(routing_matrix, dict):
            cidx = {c: i for i, c in enumerate(classes)}
            for (from_class, to_class), m in routing_matrix.items():
                arr = _as2d(m)
                if arr is not None:
                    yield cidx.get(from_class), cidx.get(to_class), arr

        elif isinstance(routing_matrix, (list, np.ndarray)):
            # P[r][s] (class switching), P[r] (per class), or a plain 2D matrix
            # applied identically to every class.
            if isinstance(routing_matrix, list) and len(routing_matrix) == nclasses \
                    and nclasses > 0 and all(isinstance(row, list) and len(row) == nclasses
                                             for row in routing_matrix):
                for r, row in enumerate(routing_matrix):
                    for s, entry in enumerate(row):
                        arr = _as2d(entry)
                        if arr is not None:
                            yield r, s, arr
                return
            if isinstance(routing_matrix, list) and len(routing_matrix) == nclasses \
                    and nclasses > 0:
                percls = [_as2d(row) for row in routing_matrix]
                if all(a is not None and a.shape == (nnodes, nnodes) for a in percls):
                    for r, arr in enumerate(percls):
                        yield r, r, arr
                    return
            arr = _as2d(routing_matrix)
            if arr is not None:
                for r in range(max(nclasses, 1)):
                    yield r, r, arr

    def _validate_routing_probabilities(self, routing_matrix) -> None:
        """
        Assert that the supplied routing probabilities are nonnegative and that
        the total leaving a node in a given class does not exceed 1.

        Both must be checked on the raw input, before conversion: every
        conversion branch in ``link`` installs an entry only ``if
        matrix_arr[i, j] > 0``, so a negative entry is silently discarded and
        the resulting model differs from the one the user wrote, with no
        diagnostic. Nonnegativity also cannot be folded into the total, because
        a negative entry can only lower a total and so passes that test
        silently; the traffic equations would then be solved over a matrix that
        is not a stochastic kernel. Matches the equivalent block in the MATLAB
        ``@MNetwork/link.m``.

        Disabled by ``model.set_checks(False)``, matching the ``enableChecks``
        gate on the equivalent MATLAB and JAR blocks. This is what exempts the
        layer models ``SolverLN`` generates, which encode an and-fork as
        simultaneous emission (a Delay feeding both Join_PreAnd and
        Fork_PostAnd at 1.0, total 2.0) and are not stochastic kernels;
        ``buildLayersRecursive`` and ``solver_ln`` both turn checks off for
        exactly that reason.

        Raises:
            ValueError: If a routing probability is negative, or if a node's
                outgoing total in some class exceeds 1.
        """
        if not getattr(self, '_do_checks', True):
            return

        # structural (non-probability) entries (Fork/SPN incidence/Router pre-strategy) skip row-sum; see _kb/04-networkstruct.md Routing-matrix validation.
        from .nodes import Fork, Place, Transition, Router

        nodes = self.get_nodes()
        classes = self.get_classes()
        nclasses = len(classes)
        nnodes = len(nodes)
        neg_tol = -GlobalConstants.FineTol

        def _name(seq, k):
            if k is None or not (0 <= k < len(seq)):
                return str(k)
            return getattr(seq[k], 'name', None) or str(k)

        blocks = list(self._iter_routing_blocks(routing_matrix))

        # (1) Nonnegativity.
        for r, s, arr in blocks:
            neg = np.argwhere(arr < neg_tol)
            if neg.size == 0:
                continue
            i, j = int(neg[0][0]), int(neg[0][1])
            if nclasses > 1 and r is not None:
                raise ValueError(
                    "[{0}] Negative routing probability {1:g} from node {2} to node {3} "
                    "(class {4} to class {5}). Routing probabilities must be "
                    "nonnegative.".format(self.name, arr[i, j], _name(nodes, i),
                                          _name(nodes, j), _name(classes, r),
                                          _name(classes, s)))
            raise ValueError(
                "[{0}] Negative routing probability {1:g} from node {2} to node {3}. "
                "Routing probabilities must be nonnegative.".format(
                    self.name, arr[i, j], _name(nodes, i), _name(nodes, j)))

        # per-node per-class outgoing total, summed over destination AND arrival class s (class switching).
        totals = {}
        for r, s, arr in blocks:
            rows = min(nnodes, arr.shape[0])
            for i in range(rows):
                totals[(i, r)] = totals.get((i, r), 0.0) + float(arr[i, :].sum())

        for (i, r), total in sorted(totals.items()):
            if total <= 1.0 + GlobalConstants.FineTol:
                continue
            if i >= nnodes:
                continue
            node = nodes[i]
            if isinstance(node, (Fork, Place, Transition, Router)):
                continue
            sched = getattr(node, 'sched_strategy', None)
            if sched is not None and sched == SchedStrategy.FORK:
                continue
            if nclasses > 1 and r is not None:
                raise ValueError(
                    "[{0}] The total routing probability for jobs leaving node {1} in "
                    "class {2} is {3:g}, which is greater than 1.0.".format(
                        self.name, _name(nodes, i), _name(classes, r), total))
            raise ValueError(
                "[{0}] The total routing probability for jobs leaving node {1} is {2:g}, "
                "which is greater than 1.0.".format(self.name, _name(nodes, i), total))

    def link(self, routing_matrix) -> None:
        """
        Set the routing matrix for the network.

        Args:
            routing_matrix: RoutingMatrix, dict with (from_class, to_class) keys,
                           or 2D list/array (applies same routing to all classes)

        Raises:
            ValueError: If routing matrix dimensions don't match network, if any
                routing probability is negative, or if a node's outgoing total
                in some class exceeds 1.
        """
        self._validate_routing_probabilities(routing_matrix)

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

            # routing_matrix format: P[r][s] (class switching), P[r] (per-class), or plain 2D (shared).

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

        self._sanitize()

        # inject deferred cache retrieval-system routing edges registered by Cache.set_retrieval_system; see _kb/09-ldes-and-cache.md.
        self._inject_retrieval_routing()

        # sync routing probabilities to node._prob_routing so routing strategy reports PROB (mirrors MATLAB link.m setProbRouting).
        self._sync_prob_routing_from_matrix()

        # Insert auto-generated ClassSwitch nodes for class switching routes
        self._insert_auto_class_switches()

        self._refresh_routing_matrix()
        self._reset_struct()

        # reducible-routing (absorbing state) check via RoutingMatrix directly, to avoid triggering refresh_struct recursion.
        if self._routing_matrix is not None:
            try:
                # Read the sparse per-(class,class) node->node routes, NOT toMatrix(): that
                # one is STATION-indexed and cached, so it goes stale as soon as link()
                # adds classes (retrieval systems, auto class switches) or non-station
                # nodes exist. Node-indexed adjacency mirrors MATLAB isRoutingErgodic.m.
                nnodes_rt = len(self._nodes)
                node_index = {id(node): i for i, node in enumerate(self._nodes)}
                adj = np.zeros((nnodes_rt, nnodes_rt))
                for routes in self._routing_matrix._routes.values():
                    for (node_src, node_dst), prob in routes.items():
                        i = node_index.get(id(node_src))
                        j = node_index.get(id(node_dst))
                        if i is not None and j is not None and prob > 0:
                            adj[i, j] = 1.0
                if nnodes_rt > 0 and np.any(adj > 0):
                    # absorbing station: routed but with no outgoing edge beyond a self-loop; unrouted stations and Sinks are exempt.
                    has_absorbing = False
                    from .nodes import Station as _Station
                    for i, node in enumerate(self._nodes):
                        if not isinstance(node, _Station):
                            continue
                        outgoing = adj[i, :].copy()
                        outgoing[i] = 0
                        if np.sum(outgoing) < 1e-10 and np.sum(adj[i, :]) > 1e-10:
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

        # Check that order-independent (OI) stations have a permutation-invariant
        # rate mu(c). Disabled by model.set_checks(False).
        if getattr(self, '_do_checks', True):
            from .nodes import Queue as _Queue
            import warnings
            oi_nodes = [nd for nd in self._nodes
                        if isinstance(nd, _Queue)
                        and getattr(nd._sched_strategy, 'name', None) == 'OI'
                        and nd._svc_rate_fun is not None]
            # only walk the class-switching graph when there is something to check
            Nvec, ntot = self._oi_check_bounds() if oi_nodes else (None, np.inf)
            for node in oi_nodes:
                # ntot caps the microstate LENGTH: with Nvec now chain-wide,
                # summing it would count each chain once per class and admit
                # microstates longer than the network can produce.
                cap = min(getattr(node, '_capacity', np.inf), ntot)
                ok, badc, partial = node.check_perm_invariance(Nvec, cap)
                if not ok:
                    raise ValueError(
                        "Order-independent (OI) station '%s' has a service rate function that is "
                        "not permutation-invariant: mu(c) differs for a reordering of the microstate "
                        "%s. Use SchedStrategy.PAS for order-dependent service, or disable this check "
                        "with model.set_checks(False)." % (node.get_name(), badc))
                # (1): per-job rates must be non-negative. Runs second because
                # permutation invariance is what makes mu a function of the
                # count vector, which is what the increment test walks.
                okm, badm, badr, partialm = node.check_rate_monotonicity(Nvec, cap)
                partial = partial or partialm
                if not okm:
                    raise ValueError(
                        "Order-independent (OI) station '%s' has a service rate function with a "
                        "negative per-job service rate: mu(c) DROPS when a class-%d job joins, at "
                        "microstate %s. An OI rate must satisfy mu(c_1..c_j) >= mu(c_1..c_{j-1}), "
                        "so that every mu_j(c) is non-negative. A processor-sharing total rate "
                        "(sum_j mu_{c_j})/n has this shape whenever classes have different rates -- "
                        "use SchedStrategy.PS for that station, or disable this check with "
                        "model.set_checks(False)." % (node.get_name(), badr, badm))
                if partial:
                    warnings.warn(
                        "[link] Order-independent (OI) station '%s': the permutation-invariance check "
                        "was only partial because the reachable population is large; a subset of "
                        "microstates was verified. To skip this check, call model.set_checks(False) "
                        "before link()." % node.get_name(), UserWarning)

    def _oi_check_bounds(self):
        """Per-class job bound, and the total, for the OI permutation-invariance check.

        The conserved quantity under class switching is the CHAIN population, not
        the per-class one: a class reaches counts up to its chain's total, and a
        class declared empty and filled only by a switch (a ClassSwitch target, a
        Cache Hit/Miss class) reaches them from a declared population of 0.
        Bounding each class by its OWN declared population therefore leaves
        reachable microstates unenumerated -- every mixed one, when the declared
        population is 0 -- and the check passes vacuously on rates that are
        plainly order-dependent.

        Union the classes over every class-switching edge and give each class its
        component's total population. Returns (Nvec, ntot), ntot being the total
        closed population, which bounds how many jobs one station can hold and so
        the microstate length. Without class switching each component is a single
        class and both fall back to what they were. An open class carries no
        population attribute, so Inf propagates through its whole chain as before.
        """
        from .nodes import Cache as _Cache, ClassSwitch as _ClassSwitch
        K = len(self._classes)
        pop = np.array([getattr(c, 'population', np.inf) for c in self._classes], dtype=float)
        index = {id(c): r for r, c in enumerate(self._classes)}

        parent = list(range(K))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            px, py = find(x), find(y)
            if px != py:
                parent[px] = py

        def union_classes(src, dst):
            r, t = index.get(id(src)), index.get(id(dst))
            if r is not None and t is not None:
                union(r, t)

        # (a) routing entries that switch class -- present when link() has not
        # rewritten them into an auto ClassSwitch node.
        if self._routing_matrix is not None:
            for key, routes in self._routing_matrix._routes.items():
                src, dst = key[0], key[1]
                if src is dst:
                    continue
                if any(prob > 0 for prob in routes.values()):
                    union_classes(src, dst)

        # (b) the switching nodes themselves, which is where link() puts it.
        for node in self._nodes:
            if isinstance(node, _ClassSwitch):
                M = getattr(node, '_switch_matrix', None)
                if M is None:
                    continue
                M = np.asarray(M, dtype=float)
                if M.ndim != 2:
                    continue
                for r in range(min(K, M.shape[0])):
                    for t in range(min(K, M.shape[1])):
                        if r != t and M[r, t] > 0:
                            union(r, t)
            elif isinstance(node, _Cache):
                for src, dst in list(getattr(node, '_hit_class', {}).items()) \
                        + list(getattr(node, '_miss_class', {}).items()):
                    union_classes(src, dst)
                for src, dsts in getattr(node, '_item_classes', {}).items():
                    for dst in dsts:
                        union_classes(src, dst)
                # keyed by (item, arriving class) -> retrieval class
                for key, dst in getattr(node, '_retrieval_classes', {}).items():
                    union_classes(key[1], dst)

        Nvec = np.empty(K)
        for r in range(K):
            root = find(r)
            Nvec[r] = np.sum(pop[[t for t in range(K) if find(t) == root]])
        return Nvec, float(np.sum(pop))

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

    def get_reducibility_info(self):
        """Reducibility structure of the routing, with a suggested repair each.

        Twin of the MATLAB ``@MNetwork/getReducibilityInfo``. Builds on
        ``is_routing_ergodic``, which is the one adjacency build and SCC
        decomposition; this adds ``isRoutingErgodic`` and ``suggestedFixes``.

        Returns:
            dict with isRoutingErgodic, isReducible, absorbingStations,
            transientStations, numSCCs and suggestedFixes.
        """
        is_erg, info = self.is_routing_ergodic()
        info = dict(info)
        info['isRoutingErgodic'] = is_erg
        info['suggestedFixes'] = []
        if is_erg:
            return info
        target = self._default_ergodic_target(info['absorbingStations'])
        for abs_name in info['absorbingStations']:
            if target is not None and target != abs_name:
                info['suggestedFixes'].append(
                    'Route jobs from %s back to %s (e.g., P{class}(%s, %s) = 1.0)'
                    % (abs_name, target, abs_name, target))
        return info

    def get_absorbing_stations(self):
        """The stations that are absorbing: once a job enters, it never leaves.

        Twin of the MATLAB ``@MNetwork/getAbsorbingStations``. The Sink is
        excluded, being legitimately absorbing in an open network.

        Returns:
            Tuple of (stations, node_indices).
        """
        _, info = self.is_routing_ergodic()
        stations = []
        idxs = []
        for name in info['absorbingStations']:
            st = self.get_station_by_name(name)
            if st is not None:
                stations.append(st)
                idxs.append(self.get_node_index(name))
        return stations, idxs

    def make_ergodic(self, target_node=None):
        """A routing matrix that makes the network ergodic.

        Twin of the MATLAB ``@MNetwork/makeErgodic``. Every absorbing station is
        redirected to ``target_node``, defaulting to the first Delay and then to
        the first non-absorbing station.

        Does NOT relink: apply the result with ``model.link(P)``. That is
        deliberate -- relinking rebuilds the struct, and the modeller should see
        the repair before it is applied.
        """
        from line_solver.api.io.logging import line_warning
        is_erg, info = self.is_routing_ergodic()
        P = self.get_linked_routing_matrix()
        if P is None:
            P = self.init_routing_matrix()
        P = [[None if blk is None else np.array(blk, dtype=float) for blk in row] for row in P]
        if is_erg:
            line_warning('makeErgodic', 'Routing is already ergodic. No changes needed.')
            return P
        if not info['absorbingStations']:
            line_warning('makeErgodic', 'No absorbing stations found, but the network is not '
                                        'ergodic. Manual intervention is required.')
            return P

        if target_node is None:
            target_name = self._default_ergodic_target(info['absorbingStations'])
        elif isinstance(target_node, str):
            target_name = target_node
        else:
            target_name = target_node.name
        if target_name is None:
            raise ValueError('makeErgodic: cannot find a suitable target node for routing.')
        # get_node_index is 1-BASED; the routing blocks are 0-based
        target_idx = self.get_node_index(target_name) - 1
        if target_idx < 0:
            raise ValueError('makeErgodic: target node "%s" not found in the network.'
                             % target_name)

        nclasses = len(self._classes)
        for abs_name in info['absorbingStations']:
            abs_idx = self.get_node_index(abs_name) - 1
            if abs_idx == target_idx:
                line_warning('makeErgodic', 'Cannot route %s to itself. Skipping.' % abs_name)
                continue
            for r in range(nclasses):
                for s in range(nclasses):
                    if P[r][s] is not None and abs_idx < P[r][s].shape[0]:
                        P[r][s][abs_idx, :] = 0.0
                if P[r][r] is not None and abs_idx < P[r][r].shape[0] \
                        and target_idx < P[r][r].shape[1]:
                    P[r][r][abs_idx, target_idx] = 1.0
        return P

    def _default_ergodic_target(self, absorbing):
        """First Delay, else the first station that is not itself absorbing."""
        from line_solver.lang.nodes import Delay as _Delay
        for node in self._nodes:
            if isinstance(node, _Delay):
                return node.name
        first = None
        for node in self._nodes:
            if node not in self._stations:
                continue
            if first is None:
                first = node.name
            if node.name not in absorbing:
                return node.name
        return first

    def _sanitize(self) -> None:
        """Preprocess the model for consistent parameterization, as MATLAB `sanitize.m`.

        Two parts, in the reference's own order.

        (1) DEFAULTS, applied to every node whose per-class parameterization is
        incomplete: an unconfigured class gets a Disabled service (Queue, Delay,
        Place, Transition) or a Disabled arrival (Source), zero capacity and the
        WAITQ drop rule, a Join admits every class without bound, a Sink routes
        nothing, and a class a station cannot serve has its routing DISABLED
        there so no visit is generated for it.

        (2) CHECKS, skipped under `set_checks(False)`: a Queue or Delay must
        serve some class, a closed class must be served somewhere, and SEPT/LEPT
        need distinct means to order by. Cache, Petri-net and fork-join models
        are exempt from the service checks because a station of theirs
        legitimately carries no per-class service.

        NOT MIRRORED, deliberately. MATLAB's Cache block materializes a default
        `accessProb` (the per-item list-transition matrix) and rounds
        `hitClass`/`missClass` to integers. Python resolves the absent access
        graph at each consumer instead (`accost is None` -> the linear chain) and
        carries hit/miss as JobClass objects rather than indices, so neither has
        state to initialize here. See _kb/11-conventions-and-gotchas.md.
        """
        from .base import DropStrategy, RoutingStrategy, SchedStrategy, Station
        from .nodes import (Cache, Delay, Fork, Join, Place, Queue, Sink, Source,
                            Transition)
        from ..distributions import Disabled

        classes = list(self._classes)

        def _disabled(dist) -> bool:
            return dist is None or isinstance(dist, Disabled)

        # ---- (1) defaults ------------------------------------------------
        for node in self._nodes:
            if isinstance(node, Join):
                for k in classes:
                    node.set_class_capacity(k, float('inf'))
                    node.set_drop_rule(k, DropStrategy.WAITQ)
            elif isinstance(node, Source):
                for k in classes:
                    if _disabled(getattr(node, '_arrival_process', {}).get(k)):
                        node.set_arrival(k, Disabled())
            elif isinstance(node, Sink):
                for k in classes:
                    node.set_routing(k, RoutingStrategy.DISABLED)
            elif isinstance(node, (Queue, Delay)):
                # Place and Transition are NOT filled here. MATLAB fills their
                # internal `server.serviceProcess` cell, which Python does not
                # model at all: a Place carries a service only once it is a
                # QUEUEING place, and a Transition's timing lives on its modes.
                for k in classes:
                    if k not in node._service_process:
                        node._service_process[k] = Disabled()
                        # the capacity is written straight to the field, as
                        # MATLAB does: set_class_capacity refuses 0, and 0 is
                        # exactly what a class with no service must be capped at
                        node._class_capacity[k] = 0
                        node.set_drop_rule(k, DropStrategy.WAITQ)
                        node._sched_param[k] = 0

        # A class a station cannot serve must not be routed to it: MATLAB does
        # this from sanitize, and it is what keeps the class out of the visit
        # ratios. Fork/Join and Petri-net stations are exempt for the same
        # reason they are exempt from the service check below.
        special = any(isinstance(nd, (Cache, Place, Transition, Fork, Join))
                      for nd in self._nodes)
        if not special:
            for node in self._nodes:
                if not isinstance(node, (Queue, Delay)):
                    continue
                if getattr(node, '_hetero_service', None):
                    continue
                for k in classes:
                    if _disabled(node._service_process.get(k)):
                        node.set_routing(k, RoutingStrategy.DISABLED)

        # Every station carries a per-class service slot, Disabled where the
        # class is not served; the Source is exempt (its slot is the arrival).
        source = next((nd for nd in self._nodes if isinstance(nd, Source)), None)
        for node in self._nodes:
            if not isinstance(node, Station) or node is source:
                continue
            if isinstance(node, Cache):
                # a Cache never re-enters itself; the read is switched, not looped
                for k in classes:
                    try:
                        node.set_prob_routing(k, node, 0.0)
                    except Exception:
                        pass
                continue
            store = getattr(node, '_service_process', None)
            if store is None:
                continue
            for k in classes:
                if k not in store:
                    store[k] = Disabled()

        if not getattr(self, '_do_checks', True):
            return

        # ---- (2) checks --------------------------------------------------
        for node in self._nodes:
            if not isinstance(node, (Queue, Delay)):
                continue
            hetero = bool(getattr(node, '_hetero_service', None))
            enabled = any(not _disabled(d) for d in node._service_process.values())
            if not enabled and not hetero and not special:
                raise ValueError("[{0}] {1} '{2}' has no service configured for any "
                                 "job class. Use setService() to configure service "
                                 "times.".format(self.name,
                                                 'Delay' if isinstance(node, Delay) else 'Queue',
                                                 node.getName()))
            self._sanitize_sept_lept(node, classes)

        if special:
            return

        # A closed class must be served somewhere; an open class may route
        # straight to the Sink, and a class read at a Cache is served there.
        stations = [nd for nd in self._nodes if isinstance(nd, Station)]
        for k in classes:
            if type(k).__name__ == 'OpenClass':
                continue
            if any(not _disabled(getattr(nd, '_service_process', {}).get(k))
                   for nd in stations):
                continue
            if any(isinstance(nd, Cache) and not _disabled(
                    getattr(nd, '_read_process', {}).get(k)) for nd in self._nodes):
                continue
            raise ValueError("[{0}] Job class '{1}' has no service configured at any "
                             "station. Every job class must have service configured at "
                             "least one station using setService().".format(
                                 self.name, k.getName()))

    def _sanitize_sept_lept(self, node, classes) -> None:
        """Rank the class service means for SEPT/LEPT, as MATLAB sanitize.m does.

        SEPT and LEPT order the classes by expected processing time, so two
        classes with the SAME mean leave the order undefined: the reference
        refuses rather than break the tie arbitrarily, and a Queue is refused
        here for the same reason. A Delay is not -- it has no queue to order,
        and the reference only sorts its parameter there.
        """
        from .nodes import Delay
        from ..distributions import Disabled

        # compared by NAME: a node carries the user-facing SchedStrategy enum,
        # which is a DIFFERENT object from lang.base's, so `is`/`==` on the
        # members silently misses. See _kb/07-cross-language-parity.md.
        sched = getattr(node, '_sched_strategy', None)
        name = getattr(sched, 'name', None)
        if name not in ('SEPT', 'LEPT'):
            return
        means = []
        for k in classes:
            d = node._service_process.get(k)
            if d is None or isinstance(d, Disabled):
                means.append(float('nan'))
            else:
                means.append(float(d.getMean()))
        arr = np.asarray(means, dtype=float)
        if isinstance(node, Delay):
            # MATLAB sorts the Delay parameter and applies no distinctness test
            for rank, idx in enumerate(np.argsort(arr, kind='stable')):
                node._sched_param[classes[idx]] = rank + 1
            return
        # unique() counts NaN once, matching MATLAB's unique on a NaN column
        finite = arr[~np.isnan(arr)]
        n_unique = len(np.unique(finite)) + (1 if np.isnan(arr).any() else 0)
        if n_unique != len(classes):
            raise ValueError('%s does not support identical service time means.' % name)
        order = np.sort(np.unique(finite))
        if name == 'LEPT':
            order = order[::-1]
        order = list(order)
        for idx, k in enumerate(classes):
            if np.isnan(arr[idx]):
                node._sched_param[k] = len(order) + 1
            else:
                node._sched_param[k] = order.index(arr[idx]) + 1

    def _inject_retrieval_routing(self) -> None:
        """
        Inject the deferred retrieval-system routing edges registered by
        Cache.set_retrieval_system into the routing matrix.

        The retrieval-pending / -complete classes are auto-generated by
        set_retrieval_system and are not part of the user-supplied routing
        matrix, so their edges are recorded on the Cache node and merged in
        here (before auto class-switch insertion processes the routing).
        """
        if self._routing_matrix is None:
            return
        from .nodes import Cache
        rm = self._routing_matrix

        def _prob_of(c1, c2, n1, n2):
            return rm._routes.get((c1, c2), {}).get((n1, n2), 0.0)

        for node in self._nodes:
            if not isinstance(node, Cache):
                continue
            qmap = getattr(node, '_retrieval_system_queue_indices', {})
            if qmap:
                cache_node = node
                n_items = node._num_items
                # (1) default per-item routing inherited from the read class topology in P
                for read_idx0, q_indices in qmap.items():
                    read_class = self._classes[int(read_idx0)]
                    q_nodes = [self._nodes[int(qi)] for qi in q_indices]
                    for it in range(n_items):
                        r_class = node.get_retrieval_class(read_class, it)
                        if r_class is None:
                            continue
                        for qn in q_nodes:
                            pe = _prob_of(read_class, read_class, cache_node, qn)  # cache -> queue
                            if pe > 0:
                                rm.set(r_class, r_class, cache_node, qn, pe)
                            px = _prob_of(read_class, read_class, qn, cache_node)  # queue -> cache
                            if px > 0:
                                rm.set(r_class, r_class, qn, cache_node, px)
                            for qn2 in q_nodes:                                    # queue -> queue
                                pqq = _prob_of(read_class, read_class, qn, qn2)
                                if pqq > 0:
                                    rm.set(r_class, r_class, qn, qn2, pqq)
                    # consume the read class's template edges over the queue set
                    for qn in q_nodes:
                        rm.set(read_class, read_class, cache_node, qn, 0.0)
                        rm.set(read_class, read_class, qn, cache_node, 0.0)
                        for qn2 in q_nodes:
                            rm.set(read_class, read_class, qn, qn2, 0.0)

            # (2) explicit deferred entries override the defaults (allow prob==0)
            for entry in getattr(node, '_retrieval_routing_entries', []):
                from_cls, to_cls, src_node, dst_node, prob = entry
                if prob < 0:
                    continue
                rm.set(
                    self._classes[int(from_cls)], self._classes[int(to_cls)],
                    self._nodes[int(src_node)], self._nodes[int(dst_node)], prob)

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

        # route through ClassSwitch: P[r,r](i,CS)=sum_s P[r,s](i,j), P[s,s](CS,j)=1.
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

                # ALL traffic on a CS edge routes through the CS node; the switch matrix alone handles the class transform.
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

        # reset routing strategies as MATLAB link.m does when class switching rebuilds dispatchers; scoped to auto-added ClassSwitch nodes only.
        self._disable_unrouted_classes()





    def _disable_unrouted_classes(self) -> None:
        """Mark DISABLED every (node, class) with no outgoing route after link.

        Mirrors the tail of MATLAB link.m, where initDispatcherJobClasses
        clears the dispatcher and setProbRouting re-adds only the classes the
        routing matrix actually sends out of the node, everything else ending
        as DISABLED. Two consumers depend on it: the struct builder must not
        hand a default RAND route to a class that cannot reach the node (a
        REPLY signal was routed into a Cache this way and rejected for having
        no item popularity), and linemodel_save serializes these strategies,
        where an absent entry is read back as RAND by the loaders.

        Strategies other than RAND/PROB (RROBIN, WRROBIN, SQ) are
        left alone: they route without appearing in the probabilistic matrix.

        The disabling is scoped to non-Station nodes, as in MATLAB, where the
        initDispatcherJobClasses loop is guarded by ~isa(nodes{ind},'Station').
        A Station keeps its default RAND dispatcher for classes absent from P:
        with class switching a job can reach a Station under a class that never
        appears as a routing source, and disabling that class strands it.
        """
        from ..constants import RoutingStrategy as _RS
        from .nodes import Sink as _Sink, Cache as _Cache, Station as _Station

        if self._routing_matrix is None:
            return
        routed = {}
        for (class_src, _class_dst), routes in self._routing_matrix._routes.items():
            for (node_src, _node_dst), prob in routes.items():
                if prob > 0:
                    routed.setdefault(id(node_src), set()).add(id(class_src))
        for node in self._nodes:
            if isinstance(node, _Sink):
                continue
            here = routed.get(id(node), set())
            for jobclass in self._classes:
                strat = getattr(node, '_routing_strategies', {}).get(jobclass)
                strat_value = None if strat is None else (strat.value if hasattr(strat, 'value') else int(strat))
                if id(jobclass) in here:
                    # a route implies an earlier DISABLED marker is stale (e.g. a REPLY signal re-enabled by LQN2QN).
                    if strat_value == _RS.DISABLED.value:
                        node.setRouting(jobclass, _RS.PROB)
                    continue
                if isinstance(node, (_Station, _Cache)):
                    # MATLAB clears dispatchers of non-station nodes only; a Cache is a station that self-switches its read class.
                    continue
                if strat_value is None or strat_value in (_RS.RAND.value, _RS.PROB.value):
                    node.setRouting(jobclass, _RS.DISABLED)

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
            self.refresh_struct()

        if self._sn is not None and hasattr(self._sn, 'rtorig') and self._sn.rtorig is not None:
            return self._sn.rtorig

        # Build rtorig from routing matrix if not available
        if self._routing_matrix is None:
            return None

        nclasses = len(self._classes)
        nnodes = len(self._nodes)

        # Create rtorig structure: list of lists of (nnodes x nnodes) matrices
        rtorig = [[np.zeros((nnodes, nnodes)) for _ in range(nclasses)] for _ in range(nclasses)]

        # _original_routes preserves pre-ClassSwitch-insertion routing, matching MATLAB sn.rtorig.
        routes_to_use = getattr(self._routing_matrix, '_original_routes', None)
        if routes_to_use is None:
            routes_to_use = self._routing_matrix._routes

        for (from_class, to_class), routes in routes_to_use.items():
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

        # logger node numbering: original nodes, then Dep-loggers, then Arv-loggers.
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

                # original routing may be a list-of-lists P[r_idx][s_idx] or a dict-of-dicts P[class_r][class_s].
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

        Checks that routing strategies are specified. Matches MATLAB
        refreshRoutingMatrix; the per-chain reference-station check lives in
        _refresh_chains, where inchain is freshly computed. Nonnegativity of
        the probabilities themselves is asserted earlier, on the raw input, by
        _validate_routing_probabilities.
        """
        if self._routing_matrix is None:
            return

        # Validate routing strategies are specified for each class
        sn = self._sn
        if sn is not None and hasattr(sn, 'routing') and sn.routing is not None:
            routing = np.asarray(sn.routing)
            K = int(sn.nclasses) if hasattr(sn, 'nclasses') else 0
            for r in range(K):
                if r < routing.shape[1] and np.all(routing[:, r] == -1):
                    import warnings
                    warnings.warn(
                        f"Routing strategy in class {r} is unspecified at all nodes.",
                        UserWarning)


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

    def get_class_names(self) -> List[str]:
        """Return the list of job-class names, indexed by class.

        References:
            MATLAB: matlab/src/lang/@MNetwork/getClassNames.m
        """
        return [cls.name for cls in self._classes]

    def get_node_names(self) -> List[str]:
        """Return the list of node names, indexed by node.

        References:
            MATLAB: matlab/src/lang/@MNetwork/getNodeNames.m
        """
        return [node.name for node in self._nodes]

    def get_station_names(self) -> List[str]:
        """Return the list of station names, indexed by station.

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (getStationNames)
        """
        return [station.name for station in self._stations]

    def get_class_switching_mask(self) -> np.ndarray:
        """Return the class-switching mask, a boolean matrix whose entry
        ``(r, s)`` is true only if jobs in class ``r`` can switch into class
        ``s`` at some node in the network.

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (getClassSwitchingMask)
        """
        return self.get_struct().csmask

    def has_product_form_solution(self) -> bool:
        """Return whether the model admits a product-form solution.

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (hasProductFormSolution)
        """
        from ..api.sn.predicates import sn_has_product_form
        return bool(sn_has_product_form(self.get_struct()))

    def get_graph(self):
        """Return a directed-graph (networkx.DiGraph) view of the network
        topology. Nodes are the network node names and edges carry the
        per-class routing probability as the ``weight`` attribute and the
        class name as the ``jobclass`` attribute.

        References:
            MATLAB: matlab/src/lang/@MNetwork/getGraph.m
        """
        import networkx as nx
        sn = self.get_struct()
        nodenames = sn.nodenames
        classnames = sn.classnames
        K = sn.nclasses
        N = sn.nnodes
        G = nx.DiGraph()
        for name in nodenames:
            G.add_node(name)
        rtnodes = sn.rtnodes
        if rtnodes is not None and getattr(rtnodes, 'size', 0) > 0:
            rtnodes = np.asarray(rtnodes)
            for i in range(N):
                for r in range(K):
                    for j in range(N):
                        for s in range(K):
                            p = rtnodes[i * K + r, j * K + s]
                            if p > 0:
                                G.add_edge(nodenames[i], nodenames[j],
                                           weight=float(p),
                                           jobclass=classnames[r] if r < len(classnames) else str(r))
        return G

    getGraph = get_graph

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

    @staticmethod
    def _dist_feature_name(dist) -> str:
        """Feature name a distribution is marked under.

        MATLAB (getUsedLangFeatures: dist.getFeatureName) and the JAR
        (Network.setUsedLangFeature(dist.getFeatureName())) both mark the name
        the distribution DECLARES, so this reads the same accessor rather than
        the Python class name. Keying on the class name meant a class whose
        declared name differed marked a feature the other codebases never mark:
        the solver gate would then fire in one codebase and not the others,
        silently, for the same model. get_feature_name is preferred over name
        because it resolves the SPECIALIZED registry entries (Cox2, Trace) that
        name cannot express; it falls back to name, then to the class name, for
        objects that define neither.
        """
        getter = getattr(dist, 'get_feature_name', None)
        if callable(getter):
            name = getter()
            if isinstance(name, str) and name:
                return name
        name = getattr(dist, 'name', None)
        if isinstance(name, str) and name:
            return name
        return type(dist).__name__

    def set_used_lang_feature(self, feature: str):
        """Set a feature as used in this model."""
        if not hasattr(self, '_used_features') or self._used_features is None:
            self.init_used_features()
        # an empty feature name means no gated capability (skipped), mirrors jar Network.getUsedLangFeatures.
        if not feature:
            return
        # unregistered process-class names now raise in set_true rather than being silently dropped from the feature gate.
        if feature in ('MarkedMAP', 'MarkedMMPP'):
            feature = 'MMAP'
        self._used_features.set_true(feature)


    def findSolver(self, metric: str = '', showAll: bool = False):
        """Which solvers and solver methods can analyze THIS model.

            model.findSolver()                # every (solver, method) pair that runs
            model.findSolver('cdf')           # ... that returns a passage-time law
            model.findSolver('getCdfRespT')   # the same question, asked by accessor
            model.findSolver('', True)        # also the pairs that are refused, and why

        The returned DataFrame has one row per pair, with columns Solver,
        Method, Runnable, Class ('exact', 'approx', 'bound' or 'simulation'),
        Metrics and Reason. Method is the method name to pass as a solver method, so
        a row can be acted on directly::

            T = model.findSolver('cdf')
            solver = LINE(model, T.Method[0])

        findMethod and help are aliases of this method.

        Args:
            metric: measure group ('cdf') or accessor ('getCdfRespT') to narrow
                the report to; '' or 'any' keeps every pair.
            showAll: also list the refused pairs, with the reason each was
                refused.

        Returns:
            pandas.DataFrame with the six columns above.
        """
        # The gate lives in SolverAUTO, which is the class that already knows
        # every family, how to build one and what each refuses. Asking it here
        # rather than reimplementing the walk is what keeps the model's answer
        # and AUTO's own dispatch from being two opinions.
        from ..solvers.solver_auto.solver_auto import SolverAUTO
        # silenced() wraps the CONSTRUCTION too: it probes every candidate
        # with supports(model), which warns on a model one of them refuses.
        with SolverAUTO.silenced():
            auto = SolverAUTO(self, verbose=False)
        return auto.findSolver(metric, showAll)

    def findMethod(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: which solvers and solver methods can analyze
        this model.

        The two names exist because the question is asked both ways round --
        "which solver do I use" and "which method do I pass" -- and the answer
        is the same table, whose Method column carries the method name either caller
        needs.
        """
        return self.findSolver(metric, showAll)

    def help(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: what can this model be solved with?"""
        return self.findSolver(metric, showAll)

    find_solver = findSolver
    find_method = findMethod

    def find_binding_capacity(self):
        """The first station whose finite capacity can actually BIND, as
        (binds, node, cap, r, is_open), or binds=False when no buffer in the
        model can refuse a job. Port of MATLAB MNetwork.findBindingCapacity.

        ONE PREDICATE, TWO CALLERS. NetworkSolver.checkBindingCapacity turns the
        answer into the refusal the product-form solvers raise, and
        get_used_lang_features marks the registry name 'FiniteCapacity' on it,
        so a solver that does not declare the name refuses exactly the models
        the structural gate refuses.

        The test reads the node-level capacity / class capacity the user set and
        the class populations from the CLASS OBJECTS, never sn.cap / sn.classcap:
        _refresh_capacity derives a FINITE classcap (the chain population) for
        every closed model, so an sn-level test would call every closed model
        capped, and reading the struct from the recorder would trigger a refresh
        on every feature query.

        Only a capacity that can bind counts. A closed model whose station
        capacity is at least the total population can never block a job, so the
        declaration is a no-op (set_capacity(N) on a station of an N-job closed
        model is a common idiom). The population of an open class is inf, so any
        finite capacity an open class can reach binds. A Cache model is exempt:
        Cache sets class_capacity=1 on the retrieval queues it builds, and the
        cache analyzers solve those rather than treating them as a buffer
        constraint.

        r is the 0-BASED class index of a per-class buffer and -1 for a
        station-level one; is_open says whether an open class reaches the
        buffer, which decides the fallback advice.
        """
        from .base import JobClassType
        nodes = getattr(self, '_nodes', []) or []
        for node in nodes:
            if type(node).__name__ == 'Cache':
                return False, None, np.inf, -1, False
        # get_number_of_jobs() returns 0 for an open class, so the population
        # vector is rebuilt here with the inf the test needs.
        njobs = np.zeros(len(self._classes))
        for i, jobclass in enumerate(self._classes):
            if jobclass.jobclass_type == JobClassType.OPEN:
                njobs[i] = np.inf
            elif jobclass.jobclass_type == JobClassType.CLOSED:
                njobs[i] = jobclass.getNumberOfJobs()
        total_jobs = float(np.sum(njobs))  # inf as soon as one class is open
        any_open = bool(np.any(np.isinf(njobs)))
        for node in nodes:
            tname = type(node).__name__
            # Mirror MATLAB's isa(node,'Station') && ~isa(node,'Source'/'Sink')
            # filter. Today only Stations carry _capacity, but do not rely on
            # that: a future node type with the attribute must not be gated here.
            is_station = any(b.__name__ == 'Station' for b in type(node).__mro__)
            if not is_station or tname in ('Source', 'Sink') or not hasattr(node, '_capacity'):
                continue
            cap = getattr(node, '_capacity', np.inf)
            if cap is not None and np.isfinite(cap) and cap >= 0 and cap < total_jobs:
                return True, node, float(cap), -1, any_open
            ccap = getattr(node, '_class_capacity', None) or {}
            for jobclass, v in ccap.items():
                if v is None or not np.isfinite(v) or v <= 0:
                    continue
                if isinstance(jobclass, int):
                    idx = jobclass
                elif hasattr(jobclass, 'get_index0'):
                    idx = jobclass.get_index0()
                else:
                    continue
                if idx is None or idx >= njobs.size:
                    continue
                if v < njobs[idx]:
                    return True, node, float(v), int(idx), bool(np.isinf(njobs[idx]))
        return False, None, np.inf, -1, False

    findBindingCapacity = find_binding_capacity

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

        # A self-looping class is a closed class whose routing returns to its own
        # reference station. It was declared by eight solvers and never marked, so
        # a solver that omitted the name (the bounds, RCAT, QNS) was never asked;
        # it is now refused unless it declares the name.
        from .classes import SelfLoopingClass as _SelfLoopingClass
        for jobclass in self._classes:
            if isinstance(jobclass, _SelfLoopingClass):
                self.set_used_lang_feature('SelfLoopingClass')
                break

        # Signal may still be an unresolved placeholder (pre-refresh_struct) so all three signal subtypes are inspected.
        from .classes import Signal, OpenSignal, ClosedSignal, SignalType, RemovalPolicy
        for jobclass in self._classes:
            if isinstance(jobclass, ClosedSignal):
                self.set_used_lang_feature('ClosedSignal')
            elif isinstance(jobclass, OpenSignal):
                self.set_used_lang_feature('OpenSignal')
            elif isinstance(jobclass, Signal):
                # unresolved placeholder: resolves by presence of a Source
                self.set_used_lang_feature(
                    'OpenSignal' if self.get_index_source_node() >= 0 else 'ClosedSignal')
            else:
                continue
            signal_type = jobclass.signal_type
            if signal_type == SignalType.NEGATIVE:
                self.set_used_lang_feature('SignalType_NEGATIVE')
            elif signal_type == SignalType.REPLY:
                self.set_used_lang_feature('SignalType_REPLY')
            elif signal_type == SignalType.CATASTROPHE:
                self.set_used_lang_feature('SignalType_CATASTROPHE')
            # removal-batch distribution and removal policy are independent capabilities, reported separately.
            if jobclass.removal_distribution is not None:
                self.set_used_lang_feature('SignalBatchRemoval')
            if (jobclass.removal_policy is not None
                    and jobclass.removal_policy != RemovalPolicy.RANDOM):
                self.set_used_lang_feature('SignalRemovalPolicy')

        # Finite capacity regions: solvers that ignore FCRs must be able to
        # reject them via the registry 'Region' feature
        if getattr(self, '_regions', None):
            self.set_used_lang_feature('Region')

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
                if getattr(node, '_service_process', None) or hasattr(node, '_service'):
                    for jobclass, dist in (getattr(node, '_service_process', None) or getattr(node, '_service', {})).items():
                        if dist is not None:
                            dist_name = Network._dist_feature_name(dist)
                            # Immediate/Disabled are internal placeholders, not
                            # user-facing distributions (same exclusion as JAR)
                            if dist_name not in ('Immediate', 'Disabled'):
                                self.set_used_lang_feature(dist_name)
                # Check routing strategy (stored on the Node base as
                # _routing_strategies: Dict[JobClass, RoutingStrategy])
                for strategy in getattr(node, '_routing_strategies', {}).values():
                    if strategy is not None:
                        self.set_used_lang_feature(RoutingStrategy.to_feature(strategy))

            elif node_type == 'Delay':
                self.set_used_lang_feature('Delay')
                self.set_used_lang_feature('SchedStrategy_INF')
                # Check service distribution
                if getattr(node, '_service_process', None) or hasattr(node, '_service'):
                    for jobclass, dist in (getattr(node, '_service_process', None) or getattr(node, '_service', {})).items():
                        if dist is not None:
                            dist_name = Network._dist_feature_name(dist)
                            # Immediate/Disabled are internal placeholders, not
                            # user-facing distributions (same exclusion as JAR)
                            if dist_name not in ('Immediate', 'Disabled'):
                                self.set_used_lang_feature(dist_name)
                # a Delay routes like any station; recorded in the same branch as Queue.
                for strategy in getattr(node, '_routing_strategies', {}).values():
                    if strategy is not None:
                        self.set_used_lang_feature(RoutingStrategy.to_feature(strategy))

            elif node_type == 'Source':
                self.set_used_lang_feature('Source')
                # Check arrival distribution
                if getattr(node, '_arrival_process', None) or hasattr(node, '_arrival'):
                    for jobclass, dist in (getattr(node, '_arrival_process', None) or getattr(node, '_arrival', {})).items():
                        if dist is not None:
                            dist_name = Network._dist_feature_name(dist)
                            # Immediate/Disabled are internal placeholders, not
                            # user-facing distributions (same exclusion as JAR)
                            if dist_name not in ('Immediate', 'Disabled'):
                                self.set_used_lang_feature(dist_name)

            elif node_type == 'Sink':
                self.set_used_lang_feature('Sink')

            elif node_type == 'Cache':
                self.set_used_lang_feature('Cache')
                self.set_used_lang_feature('CacheClassSwitcher')
                # Check replacement strategy
                if hasattr(node, '_replacement_strategy') and node._replacement_strategy is not None:
                    self.set_used_lang_feature(ReplacementStrategy.to_feature(node._replacement_strategy))
                # per-list storage cost caps with per-item sizes
                if getattr(node, '_cost_cap', None) is not None:
                    self.set_used_lang_feature('CacheItemSize')
                # delayed-hit retrieval classes are unsupported by JMT; mirrors MATLAB getUsedLangFeatures.
                if getattr(node, '_retrieval_class_indices', None):
                    self.set_used_lang_feature('CacheRetrieval')

            elif node_type == 'ClassSwitch':
                self.set_used_lang_feature('ClassSwitch')
                self.set_used_lang_feature('StatelessClassSwitcher')

            elif node_type == 'Fork':
                self.set_used_lang_feature('Fork')
                self.set_used_lang_feature('Forker')
                # variable forking levels: a name declared but never marked is a
                # name no solver can ever refuse, so mark all three here
                if getattr(node, '_tasks_per_link_by_dest', None):
                    self.set_used_lang_feature('ForkFanoutVector')
                if getattr(node, '_tasks_per_link_dist', None):
                    self.set_used_lang_feature('ForkFanoutRandom')
                if getattr(node, '_branch_prob', None):
                    self.set_used_lang_feature('ForkBranchProbability')

            elif node_type == 'Join':
                self.set_used_lang_feature('Join')
                self.set_used_lang_feature('Joiner')

            elif node_type == 'Router':
                # a Router with plain PROB/RAND is solvable by product-form solvers; only its routing strategy is recorded.
                for strategy in getattr(node, '_routing_strategies', {}).values():
                    if strategy is not None:
                        self.set_used_lang_feature(RoutingStrategy.to_feature(strategy))

            elif node_type == 'Place':
                self.set_used_lang_feature('Place')
                if getattr(node, '_queueing', False):
                    self.set_used_lang_feature('QueueingPlace')
                self.set_used_lang_feature('Storage')
                self.set_used_lang_feature('Linkage')

            elif node_type == 'Transition':
                self.set_used_lang_feature('Transition')
                self.set_used_lang_feature('Enabling')
                self.set_used_lang_feature('Timing')
                self.set_used_lang_feature('Firing')
                # Inhibitor arcs: flag only when a finite threshold is set, so
                # plain SPNs are not gated out of solvers lacking it.
                inh_conds = getattr(node, '_inhibiting_conditions', None)
                if inh_conds is not None:
                    for inh_m in inh_conds:
                        if np.any(np.isfinite(np.asarray(inh_m, dtype=float))):
                            self.set_used_lang_feature('Inhibiting')
                            break

        # A finite-server station serving several jobs at once. A Delay carries
        # inf here and is NOT one: the single-server recursions (aql, mvac, rqna,
        # explicit, the aba..scb bounds, RCAT, the M/G/1 closed forms) answered a
        # c-server station as one server of the same rate, and only a structural
        # predicate could say so.
        for node in self._nodes:
            if type(node).__name__ not in ('Queue', 'Delay'):
                continue
            nserv = getattr(node, '_number_of_servers', None)
            if nserv is not None and np.isfinite(nserv) and nserv > 1:
                self.set_used_lang_feature('MultiServer')
                break

        # Batch arrivals (Source.set_arrival_batch): a batch releases several jobs
        # at one arrival epoch, which changes the arrival stream itself. Only the
        # LDES engine reads sn.arrivalbatch; every other solver would serve the
        # single-arrival stream of the same rate, so the name is marked and left
        # to LDES to declare.
        for node in self._nodes:
            if type(node).__name__ != 'Source':
                continue
            batches = getattr(node, '_arrival_batch', None) or {}
            if any(b is not None for b in batches.values()):
                self.set_used_lang_feature('BatchArrival')
                break

        # A retrial orbit (Queue.set_retrial / set_orbit): the same per-class test
        # _refresh_balking_retrial makes for sn.retrialProc, a configured delay
        # that is not the Disabled placeholder. A solver that reads no sn.retrial*
        # field would answer the model with the refused jobs simply lost, so the
        # orbit is gated by name.
        for node in self._nodes:
            if type(node).__name__ != 'Queue':
                continue
            delays = getattr(node, '_retrial_delays', None) or {}
            if any(d is not None and type(d).__name__ != 'Disabled' for d in delays.values()):
                self.set_used_lang_feature('Retrial')
                break

        # A station or per-class buffer that can BIND: the one predicate
        # NetworkSolver.checkBindingCapacity refuses on (node-level caps against
        # the class populations, open classes always bind, Cache models exempt),
        # asked here so the refusal has a registry name and a solver method that
        # does not declare it is gated on exactly the models the structural gate
        # refuses.
        if self.find_binding_capacity()[0]:
            self.set_used_lang_feature('FiniteCapacity')

        # Limited load-dependent scaling (set_load_dependence): registered so that a
        # solver which never reads sn.lldscaling rejects the model instead of returning
        # the alpha == 1 answer. It was declared by SolverMVA/NC/CTMC/SSA/LDES and never
        # emitted here, so the gate could not fire and SolverMAM and the fluid methods
        # outside the closing family silently dropped the scaling.
        for node in self._nodes:
            get_ld = getattr(node, 'get_load_dependence', None)
            if get_ld is not None:
                lld = get_ld()
                if lld is not None and np.asarray(lld).size > 0:
                    self.set_used_lang_feature('LoadDependence')
                    break

        # class-dependence scaling registered so non-supporting solvers (e.g. JMT) reject rather than silently ignore it.
        for node in self._nodes:
            get_cd = getattr(node, 'get_class_dependence', None)
            if get_cd is not None and get_cd() is not None:
                self.set_used_lang_feature('ClassDependence')
                break

        # joint-dependence scaling (non-product-form eta_i) registered so
        # non-supporting solvers reject rather than silently ignore it.
        for node in self._nodes:
            get_jd = getattr(node, 'get_joint_dependence', None)
            if get_jd is not None and get_jd() is not None:
                self.set_used_lang_feature('JointDependence')
                break
        # Globally state-dependent scaling phi(n) over the full network state, the
        # Whittle primitive. Only SolverCTMC plumbs it, so every other solver must
        # reject the model rather than solve it unscaled.
        if self._gd_scaling is not None:
            self.set_used_lang_feature('GlobalDependence')

        # setup/delay-off registered so non-supporting solvers reject rather than silently solve setup-free.
        for node in self._nodes:
            is_enabled = getattr(node, 'is_delay_off_enabled', None)
            if is_enabled is not None and is_enabled():
                self.set_used_lang_feature('SetupDelayOff')
                break

        # A job seizing n>1 servers changes the effective capacity of the station,
        # so a solver that cannot honour it must reject the model rather than solve
        # it as if every job seized one server.
        for node in self._nodes:
            has_par = getattr(node, 'has_server_parallelism', None)
            if has_par is not None and has_par():
                self.set_used_lang_feature('ServerParallelism')
                break

        # impatience registered so non-supporting solvers (MVA/NC/MAM/FLD) reject rather than return the impatience-free rho/(1-rho).
        from .base import ImpatienceType
        for node in self._nodes:
            types = getattr(node, '_impatience_types', None)
            if types and any(t == ImpatienceType.RENEGING for t in types.values()):
                self.set_used_lang_feature('Reneging')
                break
        for node in self._nodes:
            if getattr(node, '_balking_strategies', None):
                self.set_used_lang_feature('Balking')
                break

        # breakdown registered so solvers that don't expand the (queue,server-status) chain reject rather than return the always-up result.
        for node in self._nodes:
            has_bd = getattr(node, 'has_breakdown', None)
            if has_bd is not None and has_bd():
                self.set_used_lang_feature('Breakdown')
                break

        # A genuine QUORUM join (PARTIAL, k of n siblings with k < n) fires on the
        # k-th branch completion, which is not their maximum, and the n-k stragglers
        # are discarded on arrival: a solver that serves it as a full join returns a
        # silent wrong answer, so the rule is gated by its own name. k >= n and
        # k <= 0 are full joins and are not flagged.
        from .nodes import Join as _Join, Fork as _Fork
        from .base import JoinStrategy as _JoinStrategy
        if self._nodes:
            conn = self.get_connection_matrix()
            for node_idx, node in enumerate(self._nodes):
                if not isinstance(node, _Join):
                    continue
                # siblings are counted at the FORK, as the engines do: its out-degree
                # times tasksPerLink. The join in-degree is the fallback.
                nsib = int(np.count_nonzero(conn[:, node_idx]))
                fork = node.get_fork()
                if isinstance(fork, _Fork):
                    fork_idx = self._nodes.index(fork)
                    tpl = fork.get_tasks_per_link()
                    w = 1
                    if tpl is not None:
                        tpl = np.atleast_1d(tpl)
                        if tpl.size > 0:
                            w = max(1, int(round(float(tpl.flat[0]))))
                    nsib = int(np.count_nonzero(conn[fork_idx, :])) * w
                for jobclass in self._classes:
                    if node.get_strategy(jobclass) == _JoinStrategy.STD:
                        continue
                    kreq = node.get_required(jobclass)
                    if kreq is None or kreq <= 0:
                        continue
                    if nsib <= 0 or kreq < nsib:
                        self.set_used_lang_feature('JoinPartial')

        # Heterogeneous server pools (Queue.add_server_type): the station is
        # served by several pools with their own counts, class compatibilities
        # and per-(type,class) rates. A solver that reads only sn.nservers
        # answers for a homogeneous station of the same total size, which is a
        # different system, so the pools are gated by their own name rather than
        # silently flattened.
        for node in self._nodes:
            is_hetero = getattr(node, 'is_heterogeneous', None)
            if is_hetero is not None and is_hetero():
                self.set_used_lang_feature('HeteroServers')
                break

        # A FIFO depository (Place.set_departure_discipline) releases a served
        # token only after the tokens that entered service before it, so which
        # output transitions are enabled depends on the arrival order and not
        # only on the marking. No solver implements it, hence a clean refusal
        # instead of a NORMAL depository's answer under the user's name.
        from ..constants import DepartureDiscipline as _DepartureDiscipline
        for node in self._nodes:
            disc = getattr(node, '_departure_discipline', None)
            if not disc:
                continue
            if any(d != _DepartureDiscipline.NORMAL for d in disc.values()):
                self.set_used_lang_feature('DepartureDiscipline')
                break

        return self._used_features

    # MATLAB-compatible aliases
    initUsedFeatures = init_used_features
    setUsedLangFeature = set_used_lang_feature
    getUsedLangFeatures = get_used_lang_features

    # =====================================================================
    # STRUCT COMPILATION METHODS
    # =====================================================================

    def _resolve_signals(self) -> None:
        """Resolve Signal placeholders to their concrete open or closed form.

        Mirrors MATLAB @MNetwork/resolveSignals.m, which runs at the head of
        refreshStruct. A Signal is declared OPEN by default and only becomes a
        closed class once the network is known to have no Source, so without
        this step a closed model carrying a REPLY signal fails validation with
        "Open classes require Source node".

        The placeholder is resolved in place rather than replaced by a
        ClosedSignal instance: every per-class dictionary on the nodes
        (service, routing, class switching) is keyed by the class object, so
        substituting a new object would silently orphan all of them.
        """
        from .classes import Signal, OpenSignal, ClosedSignal, JobClassType

        pending = [c for c in self._classes
                   if isinstance(c, Signal)
                   and not isinstance(c, (OpenSignal, ClosedSignal))
                   and not getattr(c, '_signal_resolved', False)]
        if not pending:
            return

        is_open = self.get_index_source_node() >= 0
        if is_open:
            for sig in pending:
                sig._signal_resolved = True
            return

        default_refstat = None
        for station in self._stations:
            if station.__class__.__name__ == 'Delay':
                default_refstat = station
                break
        if default_refstat is None:
            for station in self._stations:
                if station.__class__.__name__ == 'Queue':
                    default_refstat = station
                    break
        if default_refstat is None and self._stations:
            default_refstat = self._stations[0]

        for sig in pending:
            target = getattr(sig, '_target_job_class', None)
            refstat = getattr(target, '_refstat', None) if target is not None else None
            if refstat is None:
                refstat = default_refstat
            sig._jobclass_type = JobClassType.CLOSED
            sig._njobs = 0
            sig._refstat = refstat
            sig._signal_resolved = True
            sig._invalidate_java()

    def refresh_struct(self) -> None:
        """
        Compile model to NetworkStruct for solver consumption.

        Always performs a full rebuild of the struct from scratch.
        """
        self._resolve_signals()
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
            self._rates_dirty = False
            return
        self._refresh_rates()
        self._rates_dirty = False

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
        if not self._has_struct or self._sn is None:
            self.refresh_struct()
        elif self._rates_dirty:
            # set_service marked the rates stale but kept the struct, since the
            # swap was structurally neutral. Bring them up to date lazily.
            self.refresh_rates()

        return self._sn

    def getStateSpace(self, *args):
        """Get the CTMC state space of the model.

        Delegates to SolverCTMC; results are not cached. For repeated use,
        build a SolverCTMC(model, ...) and call getStateSpace(...) on it.

        Returns:
            (stateSpace, nodeStateSpace)
        """
        from ..solvers.solver_ctmc.solver_ctmc import SolverCTMC
        return SolverCTMC(self).getStateSpace(*args)

    def get_state_space(self, *args):
        """Get the CTMC state space (snake_case alias for getStateSpace)."""
        return self.getStateSpace(*args)

    def stateSpace(self, *args):
        """Get the CTMC state space (alias for getStateSpace)."""
        return self.getStateSpace(*args)

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

        _console.compiling(str(self.getName()))
        _console.compile_detail('reading the routing strategies of %s and %s',
                                _console._plural(self._sn.nnodes, 'node'),
                                _console._plural(self._sn.nclasses, 'class', 'classes'))

        # Build node mappings
        self._refresh_node_mappings()

        # Extract service/arrival rates
        _console.compile_detail('refreshing service and arrival processes')
        self._refresh_rates()

        # Extract load-dependent scaling
        self._refresh_load_dependence()

        # Extract class-dependent scaling
        self._refresh_class_dependence()

        # Extract scheduling strategies
        _console.compile_detail('refreshing the scheduling strategies')
        self._refresh_scheduling()

        # Build routing matrix
        _console.compile_detail('computing the routing table')
        self._refresh_routing()

        # Compute chains and visits
        _console.compile_detail('computing the chains and the visit ratios')
        self._refresh_chains()
        _console.compile_detail('found %s over %s',
                                _console._plural(self._sn.nchains, 'chain'),
                                _console._plural(self._sn.nclasses, 'class', 'classes'))

        # Extract capacity (must be after _refresh_chains since capacity depends on chain info)
        self._refresh_capacity()

        # Extract node parameters (for transitions, caches, etc.)
        _console.compile_detail('refreshing node parameters and state-dependent routing')
        self._refresh_nodeparam()

        # state-dependent routing params exposed in sn.nodeparam; must run after _refresh_nodeparam rebuilds it.
        self._refresh_statedep_routing_params()
        self._refresh_state_dep_routing()

        # Build fork-join relationship matrix
        self._refresh_fork_joins()

        # struct marked built before fork-join nodevisits adjustment, to break mmt() -> get_linked_routing_matrix() -> _build_struct_from_scratch() recursion.
        self._has_struct = True

        # MMT nodevisits adjustment skipped on FJ tag-augmented copies (ModelAdapter.fjtag already overwrote them).
        if not getattr(self, 'is_fj_augmented', False):
            self._refresh_fork_join_nodevisits()

        # Extract initial state (for SPNs)
        self._refresh_state()

        # Populate finite capacity region information
        _console.compile_detail('refreshing finite capacity regions')
        self._refresh_regions()

        # Populate balking, retrial, and varsparam fields
        self._refresh_balking_retrial()
        self._refresh_varsparam()

        # heterogeneous-server fields populated last so _refresh_nodeparam cannot wipe them; mirrors MATLAB/JAR ordering.
        self._refresh_heterogeneous_servers()

        # Synchronous calls (REPLY signals): sn.syncreply and the per-node
        # blocked-server declaration sn.replyblock. Needs rtnodes and sched, so
        # it runs after routing and scheduling.
        self._refresh_syncreply()

        # Local-variable widths. Runs after scheduling and nodeparam, both of
        # which the polling controller's width is read from.
        self._refresh_local_vars()

        # Needs both classprio and sched, so it runs last.
        self._warn_priorities_ignored()

        # Needs the visit ratios, so it runs after _refresh_chains.
        self._check_service_reachable()

    #: The policies that read sn.classprio. Priority-awareness is a property of
    #: the DECLARED policy, never inferred from the data: a base policy is never
    #: upgraded because the classes it was handed carry unequal priorities.
    #: See _kb/11-conventions-and-gotchas.md.
    _PRIO_SCHEDS = (
        SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
        SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO,
        SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO,
        SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO,
        SchedStrategy.SRPTPRIO,
    )

    def _warn_priorities_ignored(self) -> None:
        """Warn when class priorities are declared but no station honors them.

        Twin of the block at the tail of MATLAB @MNetwork/refreshStruct.m and of
        Network.refreshStruct in the JAR: unequal priorities with no *PRIO
        station mean the priorities change nothing, which is worth saying once
        rather than letting the user read the resulting metrics as priority ones.

        References:
            MATLAB: matlab/src/lang/@MNetwork/refreshStruct.m
        """
        from ..api.io.logging import line_printf, line_warning

        classprio = getattr(self._sn, 'classprio', None)
        if classprio is None or len(classprio) == 0:
            return
        classprio = np.asarray(classprio).flatten()
        if np.all(classprio == classprio[0]):
            return

        sched = getattr(self._sn, 'sched', None) or {}
        prio_values = set(s.value for s in Network._PRIO_SCHEDS)
        declared = set(s.value if hasattr(s, 'value') else s for s in sched.values())
        if not (declared & prio_values):
            line_warning('refresh_struct',
                         'Priority classes are specified but no priority-aware scheduling policy '
                         '(PSPRIO, DPSPRIO, GPSPRIO, HOL, FCFSPRIO, FCFSPRPRIO, FCFSPIPRIO, '
                         'LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO, SRPTPRIO) is used in the model. '
                         'Priorities will be ignored.')
            return

        # A LOWER classprio value is the more urgent one in LINE.
        names = getattr(self._sn, 'classnames', None) or \
            ['Class%d' % (j + 1) for j in range(len(classprio))]
        high = ','.join(names[j] for j in np.flatnonzero(classprio == classprio.min()))
        low = ','.join(names[j] for j in np.flatnonzero(classprio == classprio.max()))
        line_printf('Priority: highest=%s, lowest=%s\n', high, low)

    def _refresh_syncreply(self) -> None:
        """Populate sn.syncreply and sn.replyblock (synchronous calls).

        sn.syncreply[r] is the 0-based index of the REPLY signal class that
        class r waits for, -1 when class r makes no synchronous call (the
        bidirectional link is established by Signal.forJobClass, and only for
        SignalType.REPLY).

        sn.replyblock[ind, r] is 1 when node ind holds a server for a class-r
        job that has left for its callee and is waiting for the reply. The
        holding station is identified STRUCTURALLY, since a CTMC has no job
        identity to key on as LDES does: it is a station that the REPLY class
        is routed INTO. In the canonical shape Client -> Server -> (switch to
        Reply) -> Client, that selects the client and NOT the server -- the
        server's departure is the one that CREATES the reply, and LDES likewise
        does not block it (Solver_ssj: !classSwitchedToReply). Stations that
        never receive the reply class carry no counter and are untouched.

        Mirrors MATLAB @MNetwork/refreshStruct.m (syncreply) and
        @MNetwork/refreshLocalVars.m (replyblock).
        """
        sn = self._sn
        R = int(sn.nclasses)
        nnodes = int(sn.nnodes)

        syncreply = np.full(R, -1, dtype=int)
        for r, jobclass in enumerate(self._classes):
            reply = None
            if hasattr(jobclass, 'get_reply_signal_class'):
                reply = jobclass.get_reply_signal_class()
            elif hasattr(jobclass, '_reply_signal_class'):
                reply = jobclass._reply_signal_class
            if reply is None:
                continue
            for s, cc in enumerate(self._classes):
                if cc is reply:
                    syncreply[r] = s
                    break
        sn.syncreply = syncreply

        replyblock = np.zeros((nnodes, R), dtype=int)
        sn.replyblock = replyblock
        if not np.any(syncreply >= 0):
            return

        rtnodes = np.asarray(sn.rtnodes) if sn.rtnodes is not None else None
        if rtnodes is None or rtnodes.size == 0:
            return

        def _enum(x):
            return int(x.value) if hasattr(x, 'value') else int(x)

        for r in range(R):
            s = int(syncreply[r])
            if s < 0 or s >= R:
                continue
            for ind in range(nnodes):
                if not bool(sn.isstation[ind]):
                    continue
                if _enum(sn.nodetype[ind]) == _enum(NodeType.SOURCE):
                    continue
                ist = int(sn.nodeToStation[ind])
                # An INF station has a server for every job, so holding one is
                # immaterial and needs no state.
                if _enum(sn.sched[ist]) == _enum(SchedStrategy.INF):
                    continue
                # Does the reply class ever arrive here?
                if not np.any(rtnodes[:, ind * R + s] > 0):
                    continue
                if _enum(sn.sched[ist]) != _enum(SchedStrategy.FCFS):
                    raise RuntimeError(
                        "Synchronous calls (REPLY signals) are supported only at FCFS stations, "
                        "but %s uses %s. A held server is encoded as a per-class counter, which is "
                        "exact only where servers are interchangeable (FCFS) or unlimited (INF); the "
                        "other disciplines are not yet encoded rather than infeasible. Set this station "
                        "to FCFS or INF, or simulate the layered model directly with SolverLDES."
                        % (sn.nodenames[ind], SchedStrategy(sn.sched[ist]).name.lower()))
                replyblock[ind, r] = 1

    def _refresh_state_dep_routing(self) -> None:
        """Build sn.sdr, the station-indexed Krzesinski SDR structure.

        The node-indexed twin stays in sn.nodeparam[ind][r]['sdr'], where the
        routing closure of refresh_sync reads it. A network admits one
        subnetwork Q(V,V): the routing probabilities of Krzesinski (1987) are
        chain independent, so every class routed by it must declare the same
        branches, nesting and coefficients.
        """
        import numpy as _np
        from .base import RoutingStrategy
        sdr_val = int(RoutingStrategy.SDR)
        self._sn.sdr = None
        decl, declnode = None, None
        for i, node in enumerate(self._nodes):
            strategies = getattr(node, '_routing_strategies', None)
            if not strategies:
                continue
            for jobclass, strategy in strategies.items():
                sv = strategy.value if hasattr(strategy, 'value') else int(strategy)
                if sv != sdr_val:
                    continue
                params = getattr(node, '_routing_params', {}).get(jobclass, ())
                if not params or not isinstance(params[0], dict):
                    raise ValueError('node %s declares SDR routing without a '
                                     'structure; use set_state_dep_routing'
                                     % self._sn.nodenames[i])
                cand = dict(params[0])
                cand['entry'] = node
                if decl is None:
                    decl, declnode = cand, i
                elif not _sdr_same(decl, cand):
                    raise ValueError(
                        'two different state-dependent routing structures are declared '
                        '(nodes %s and %s); the routing probabilities of Krzesinski (1987) '
                        'are chain independent, so a network admits one subnetwork Q(V,V)'
                        % (self._sn.nodenames[declnode], self._sn.nodenames[i]))
                class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                nodeidx = {
                    'entry': i,
                    'departure': self._nodes.index(cand['departure']),
                    'branch': [None] + [[self._nodes.index(x) for x in cand['branch'][b]]
                                        for b in range(1, len(cand['branch']))],
                    'entryOf': [0] + [self._nodes.index(cand['entryOf'][b])
                                      for b in range(1, len(cand['branch']))],
                    'departureOf': [0] + [self._nodes.index(cand['departureOf'][b])
                                          for b in range(1, len(cand['branch']))],
                    'level': list(cand['level']),
                    'C': list(cand['C']),
                    'd': _np.asarray(cand['d'], dtype=float),
                }
                if self._sn.nodeparam is None:
                    self._sn.nodeparam = {}
                if i not in self._sn.nodeparam or self._sn.nodeparam[i] is None:
                    self._sn.nodeparam[i] = {}
                entry = self._sn.nodeparam[i]
                if isinstance(entry, dict):
                    if not isinstance(entry.get(class_idx), dict):
                        entry[class_idx] = {}
                    entry[class_idx]['sdr'] = nodeidx
        if decl is None:
            return

        def nd2station(ind):
            ist = int(self._sn.nodeToStation[ind])
            if ist < 0:
                raise ValueError('node %s takes part in state-dependent routing but '
                                 'is not a station: the product form is over queue '
                                 'lengths, and a stateless node holds none'
                                 % self._sn.nodenames[ind])
            return ist

        B = len(decl['branch'])
        sdr = {
            'entry': nd2station(self._nodes.index(decl['entry'])),
            'departure': nd2station(self._nodes.index(decl['departure'])),
            'branch': [None] + [[nd2station(self._nodes.index(x)) for x in decl['branch'][b]]
                                for b in range(1, B)],
            'entryOf': [0] + [nd2station(self._nodes.index(decl['entryOf'][b])) for b in range(1, B)],
            'departureOf': [0] + [nd2station(self._nodes.index(decl['departureOf'][b])) for b in range(1, B)],
            'level': list(decl['level']),
            'C': list(decl['C']),
            'd': _np.asarray(decl['d'], dtype=float),
        }
        from ..api.pfqn.sdr import pfqn_sdrcoeff
        pfqn_sdrcoeff(sdr)  # validates the declaration and its population bounds
        self._sn.sdr = sdr

    def _refresh_statedep_routing_params(self) -> None:
        """Expose per-class SQ routing parameters in sn.nodeparam.

        Mirrors MATLAB, where refreshRoutingMatrix closures read d/m
        from sn.nodeparam{ind}{r}. The python refresh_sync closures read the
        same data from sn.nodeparam[ind][r]:
          SQ: 'd'
        """
        from .base import RoutingStrategy
        kch_val = int(RoutingStrategy.SQ)
        wrr_val = int(RoutingStrategy.WRROBIN)
        for i, node in enumerate(self._nodes):
            strategies = getattr(node, '_routing_strategies', None)
            if not strategies:
                continue
            for jobclass, strategy in strategies.items():
                sv = strategy.value if hasattr(strategy, 'value') else int(strategy)
                if sv not in (kch_val, wrr_val):
                    continue
                class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                params = getattr(node, '_routing_params', {}).get(jobclass, ())
                if not isinstance(params, tuple):
                    params = (params,)
                if self._sn.nodeparam is None:
                    self._sn.nodeparam = {}
                if i not in self._sn.nodeparam or self._sn.nodeparam[i] is None:
                    self._sn.nodeparam[i] = {}
                entry = self._sn.nodeparam[i]
                if not isinstance(entry, dict):
                    continue  # node already carries a structured param (e.g. Transition)
                if not isinstance(entry.get(class_idx), dict):
                    entry[class_idx] = {}
                if sv == wrr_val:
                    # WRROBIN cyclic destination list: distinct destinations repeated by integer weight; mirrors MATLAB refreshRoutingMatrix nodeparam.
                    weights_map = getattr(node, '_routing_weights', {}).get(jobclass, {})
                    idx_weight = {}
                    for dest_obj, wt in weights_map.items():
                        di = self._nodes.index(dest_obj) if dest_obj in self._nodes else int(dest_obj)
                        idx_weight[di] = int(round(float(wt))) if wt is not None else 1
                    outlinks = sorted(idx_weight.keys())
                    weighted = []
                    for di in outlinks:
                        weighted.extend([di] * max(1, idx_weight[di]))
                    entry[class_idx]['outlinks'] = outlinks
                    entry[class_idx]['weighted_outlinks'] = weighted
                else:
                    if len(params) >= 1 and params[0] is not None:
                        entry[class_idx]['d'] = int(params[0])

    def _refresh_heterogeneous_servers(self) -> None:
        """Populate per-station heterogeneous-server fields into sn.nodeparam.

        Mirrors MATLAB @MNetwork/refreshStruct.m: for each Queue with server
        types, store nservertypes, servertypenames, serverspertype,
        servercompat (nTypes x nclasses), and heteroschedpolicy on the node's
        nodeparam entry. Must run after _refresh_nodeparam.
        """
        from .nodes import Queue
        sn = self._sn
        if sn is None or sn.stationToNode is None:
            return
        nclasses = int(sn.nclasses)
        # Job parallelism is independent of the pools: a homogeneous station may
        # declare it, and a heterogeneous one may not.
        for ist, station in enumerate(self._stations):
            if not isinstance(station, Queue) or not station.has_server_parallelism():
                continue
            if sn.nodeparam is None:
                sn.nodeparam = {}
            node_idx = int(sn.stationToNode[ist])
            param = sn.nodeparam.get(node_idx)
            if param is None:
                param = {}
                sn.nodeparam[node_idx] = param
            par = np.ones(nclasses)
            for r in range(nclasses):
                par[r] = max(1, station.get_server_parallelism(self._classes[r]))
            if isinstance(param, dict):
                param['serverparallelism'] = par
            else:
                setattr(param, 'serverparallelism', par)
        for ist, station in enumerate(self._stations):
            if not isinstance(station, Queue) or not station.is_heterogeneous():
                continue
            if sn.nodeparam is None:
                sn.nodeparam = {}
            server_types = station.get_server_types()
            n_types = len(server_types)
            if n_types == 0:
                continue
            node_idx = int(sn.stationToNode[ist])
            param = sn.nodeparam.get(node_idx)
            if param is None or not isinstance(param, dict):
                if param is None:
                    param = {}
                    sn.nodeparam[node_idx] = param
            compat = np.zeros((n_types, nclasses))
            names = []
            spt = np.zeros(n_types)
            # heterorates(t,r) = service rate of class r on a type-t server
            # (1/mean; 0 where incompatible or unset - falls back to base rate).
            heterorates = np.zeros((n_types, nclasses))
            for t, st in enumerate(server_types):
                names.append(st.get_name())
                spt[t] = st.get_num_of_servers()
                for r in range(nclasses):
                    if st.is_compatible(self._classes[r]):
                        compat[t, r] = 1.0
                        dist = station.get_hetero_service(self._classes[r], st)
                        if dist is not None:
                            mval = dist.get_mean()
                            if mval > 0:
                                heterorates[t, r] = 1.0 / mval
            policy = station.get_hetero_sched_policy()
            fields = {
                'nservertypes': n_types,
                'servertypenames': names,
                'serverspertype': spt,
                'servercompat': compat,
                'heterorates': heterorates,
            }
            if policy is not None:
                fields['heteroschedpolicy'] = policy
            if isinstance(param, dict):
                param.update(fields)
            else:
                for k, v in fields.items():
                    setattr(param, k, v)

    def _refresh_balking_retrial(self) -> None:
        """Populate sn.balkingStrategy/Thresholds and sn.retrialType/Mu/Phi/Proc/MaxAttempts
        from node objects. Mirrors MATLAB refreshStruct.m (post-impatience block) and
        JAR Network.java refreshBalking()/refreshRetrial()."""
        sn = self._sn
        M = sn.nstations
        K = sn.nclasses

        sn.impatienceClass = np.zeros((M, K), dtype=int)
        sn.impatienceType = np.zeros((M, K), dtype=int)
        sn.impatienceMu = np.zeros((M, K), dtype=float)
        sn.impatiencePhi = np.zeros((M, K), dtype=float)
        sn.impatiencePhases = np.zeros((M, K), dtype=int)
        sn.impatienceProc = [[None for _ in range(K)] for _ in range(M)]
        sn.impatiencePie = [[None for _ in range(K)] for _ in range(M)]
        sn.impatienceDist = [[None for _ in range(K)] for _ in range(M)]
        sn.balkingStrategy = np.zeros((M, K), dtype=int)
        sn.balkingThresholds = [[None for _ in range(K)] for _ in range(M)]
        sn.retrialType = np.zeros((M, K), dtype=int)
        sn.retrialMu = np.zeros((M, K), dtype=float)
        sn.retrialPhi = np.zeros((M, K), dtype=float)
        sn.retrialProc = [[None for _ in range(K)] for _ in range(M)]
        sn.retrialMaxAttempts = -np.ones((M, K), dtype=int)
        from .base import RetrialPolicy as _RetrialPolicy
        sn.retrialPolicy = int(_RetrialPolicy.LINEAR) * np.ones((M, K), dtype=int)
        sn.orbitMaxJobs = -np.ones((M, K), dtype=int)
        sn.orbitImpatience = [[None for _ in range(K)] for _ in range(M)]
        sn.batchRejectProb = np.zeros((M, K), dtype=float)
        sn.impatienceType = np.zeros((M, K), dtype=int)
        sn.impatienceMu = np.zeros((M, K), dtype=float)

        # breakdown/repair are per-station (server property); degraded service is per (station,class).
        sn.hasbreakdown = np.zeros(sn.nnodes, dtype=int)
        sn.breakdownMu = np.zeros(M, dtype=float)
        sn.repairMu = np.zeros(M, dtype=float)
        sn.breakdownProc = [None for _ in range(M)]
        sn.repairProc = [None for _ in range(M)]
        sn.downServiceRates = np.zeros((M, K), dtype=float)
        self._refresh_breakdown(sn)

        _proc_type_map = {
            'Exp': ProcessType.EXP, 'Erlang': ProcessType.ERLANG,
            'HyperExp': ProcessType.HYPEREXP, 'Det': ProcessType.DET,
            'Gamma': ProcessType.GAMMA, 'Pareto': ProcessType.PARETO,
            'Weibull': ProcessType.WEIBULL, 'Lognormal': ProcessType.LOGNORMAL,
            'Uniform': ProcessType.UNIFORM, 'Coxian': ProcessType.COXIAN,
            'APH': ProcessType.APH, 'PH': ProcessType.PH,
        }

        from .nodes import Queue
        for i, station in enumerate(self._stations):
            if not isinstance(station, Queue):
                continue
            for r, jc in enumerate(self._classes):
                # impatience class set independently of whether a patience distribution is configured; mirrors MATLAB refreshStruct.
                itype = None
                if hasattr(station, '_impatience_types'):
                    itype = station._impatience_types.get(jc)
                if itype is not None:
                    sn.impatienceClass[i, r] = int(itype.value if hasattr(itype, 'value') else itype)

                # Impatience (reneging/patience)
                pdist = None
                if hasattr(station, '_patience_distributions'):
                    pdist = station._patience_distributions.get(jc)
                if pdist is not None and type(pdist).__name__ != 'Disabled':
                    cname_p = type(pdist).__name__
                    pt_p = {
                        'Exp': ProcessType.EXP, 'Erlang': ProcessType.ERLANG,
                        'HyperExp': ProcessType.HYPEREXP, 'Det': ProcessType.DET,
                        'Gamma': ProcessType.GAMMA, 'Pareto': ProcessType.PARETO,
                        'Weibull': ProcessType.WEIBULL, 'Lognormal': ProcessType.LOGNORMAL,
                        'Uniform': ProcessType.UNIFORM, 'Coxian': ProcessType.COXIAN,
                        'APH': ProcessType.APH, 'PH': ProcessType.PH,
                    }.get(cname_p, ProcessType.PH)
                    sn.impatienceType[i, r] = int(pt_p.value if hasattr(pt_p, 'value') else pt_p)
                    get_mean = getattr(pdist, 'getMean', None) or getattr(pdist, 'get_mean', None)
                    if get_mean is not None:
                        m_p = get_mean()
                        sn.impatienceMu[i, r] = (1.0 / m_p) if m_p and m_p > 0 else float('inf')
                    get_scv = getattr(pdist, 'getSCV', None) or getattr(pdist, 'get_scv', None)
                    sn.impatiencePhi[i, r] = get_scv() if get_scv is not None else 1.0
                    get_phases = getattr(pdist, 'getNumberOfPhases', None)
                    if get_phases is not None:
                        try:
                            sn.impatiencePhases[i, r] = int(get_phases())
                        except NotImplementedError:
                            sn.impatiencePhases[i, r] = 1
                    get_d0 = getattr(pdist, 'getD0', None)
                    get_d1 = getattr(pdist, 'getD1', None)
                    if get_d0 is not None and get_d1 is not None:
                        try:
                            sn.impatienceProc[i][r] = (np.array(get_d0()), np.array(get_d1()))
                        except NotImplementedError:
                            pass
                    get_pie = getattr(pdist, 'getInitProb', None)
                    if get_pie is not None:
                        try:
                            sn.impatiencePie[i][r] = np.array(get_pie())
                        except NotImplementedError:
                            pass
                    sn.impatienceDist[i][r] = pdist
                # Balking
                if hasattr(station, 'has_balking') and station.has_balking(jc):
                    strategy, thresholds = station.get_balking(jc)
                    if strategy is not None:
                        sval = strategy.value if hasattr(strategy, 'value') else strategy
                        sn.balkingStrategy[i, r] = int(sval)
                    sn.balkingThresholds[i][r] = thresholds
                # orbit shape defaults to the classical unbounded per-customer orbit unless set_orbit overrides it.
                if hasattr(station, '_retrial_policies'):
                    pol = station._retrial_policies.get(jc)
                    if pol is not None and int(pol) > 0:
                        sn.retrialPolicy[i, r] = int(pol)
                    omax = station._orbit_max_jobs.get(jc)
                    if omax is not None:
                        sn.orbitMaxJobs[i, r] = int(omax)
                # Retrial
                rdist = None
                rmax = -1
                if hasattr(station, '_retrial_delays'):
                    rdist = station._retrial_delays.get(jc)
                    if hasattr(station, '_retrial_max_attempts'):
                        rmax = station._retrial_max_attempts.get(jc, -1)
                if rdist is not None:
                    # Determine ProcessType from the distribution class
                    cname = type(rdist).__name__
                    type_map = {
                        'Exp': ProcessType.EXP, 'Erlang': ProcessType.ERLANG,
                        'HyperExp': ProcessType.HYPEREXP, 'Det': ProcessType.DET,
                        'Gamma': ProcessType.GAMMA, 'Pareto': ProcessType.PARETO,
                        'Weibull': ProcessType.WEIBULL, 'Lognormal': ProcessType.LOGNORMAL,
                        'Uniform': ProcessType.UNIFORM, 'Coxian': ProcessType.COXIAN,
                        'APH': ProcessType.APH, 'PH': ProcessType.PH,
                    }
                    pt = type_map.get(cname, ProcessType.PH)
                    sn.retrialType[i, r] = int(pt.value if hasattr(pt, 'value') else pt)
                    get_mean = getattr(rdist, 'getMean', None) or getattr(rdist, 'get_mean', None)
                    if get_mean is not None:
                        m = get_mean()
                        sn.retrialMu[i, r] = (1.0 / m) if m and m > 0 else float('inf')
                    get_scv = getattr(rdist, 'getSCV', None) or getattr(rdist, 'get_scv', None)
                    sn.retrialPhi[i, r] = get_scv() if get_scv is not None else 1.0
                    get_proc = getattr(rdist, 'getProcess', None) or getattr(rdist, 'get_process', None)
                    if get_proc is not None:
                        sn.retrialProc[i][r] = get_proc()
                    else:
                        # Native distributions may lack getProcess but expose
                        # getD0/getD1 returning the renewal MAP representation
                        get_d0 = getattr(rdist, 'getD0', None)
                        get_d1 = getattr(rdist, 'getD1', None)
                        if get_d0 is not None and get_d1 is not None:
                            sn.retrialProc[i][r] = (np.array(get_d0()), np.array(get_d1()))
                    sn.retrialMaxAttempts[i, r] = int(rmax)
                # Orbit impatience (abandonment from the retrial orbit); store (D0,D1) process
                odist = None
                if hasattr(station, '_orbit_impatience'):
                    odist = station._orbit_impatience.get(jc)
                if odist is not None and type(odist).__name__ != 'Disabled':
                    get_d0 = getattr(odist, 'getD0', None)
                    get_d1 = getattr(odist, 'getD1', None)
                    if get_d0 is not None and get_d1 is not None:
                        sn.orbitImpatience[i][r] = (np.array(get_d0()), np.array(get_d1()))
                # Batch rejection probability (retrial queues: probability that an
                # arriving batch is rejected in full when it does not fit)
                if hasattr(station, '_batch_reject_prob'):
                    sn.batchRejectProb[i, r] = float(station._batch_reject_prob.get(jc, 0.0))
                # Reneging (queue patience): impatienceType = ProcessType of the
                # patience distribution, impatienceMu = 1/mean. Only RENEGING type.
                pdist = None
                if hasattr(station, '_patience_distributions'):
                    pdist = station._patience_distributions.get(jc)
                if pdist is not None and type(pdist).__name__ != 'Disabled':
                    from .base import ImpatienceType
                    itype = (station._impatience_types.get(jc)
                             if hasattr(station, '_impatience_types') else ImpatienceType.RENEGING)
                    if itype == ImpatienceType.RENEGING:
                        pt = _proc_type_map.get(type(pdist).__name__, ProcessType.PH)
                        sn.impatienceType[i, r] = int(pt.value if hasattr(pt, 'value') else pt)
                        get_mean = getattr(pdist, 'getMean', None) or getattr(pdist, 'get_mean', None)
                        if get_mean is not None:
                            m = get_mean()
                            sn.impatienceMu[i, r] = (1.0 / m) if m and m > 0 else float('inf')

    def _refresh_breakdown(self, sn) -> None:
        """Populate the server breakdown / repair fields of sn.

        Mirrors the breakdown block of MATLAB @MNetwork/refreshStruct.m. Failure
        and repair are properties of the SERVER, so they are per station; the
        optional degraded service used while the server is down is a service
        distribution and is therefore per class.
        """
        from .nodes import Queue
        from ..distributions.continuous import Exp
        from ..distributions.base import Distribution
        K = sn.nclasses
        for ist, node in enumerate(self._stations):
            if not isinstance(node, Queue) or not node.has_breakdown():
                continue
            fdist, rdist, dsvc_list = node.get_breakdown()
            sn.hasbreakdown[int(sn.stationToNode[ist])] = 1
            sn.breakdownMu[ist] = 1.0 / fdist.getMean()
            sn.repairMu[ist] = 1.0 / rdist.getMean()
            sn.breakdownProc[ist] = (np.array(fdist.getD0()), np.array(fdist.getD1()))
            sn.repairProc[ist] = (np.array(rdist.getD0()), np.array(rdist.getD1()))
            for r in range(K):
                dsvc = None
                if len(dsvc_list) == 1:
                    dsvc = dsvc_list[0]
                elif len(dsvc_list) > r:
                    dsvc = dsvc_list[r]
                if dsvc is None or not isinstance(dsvc, Distribution) or \
                        type(dsvc).__name__ == 'Disabled':
                    continue
                if not isinstance(dsvc, Exp):
                    raise ValueError(
                        "Station '%s': the down-server service distribution must be "
                        "exponential; a phase-type degraded service would need its own "
                        "phase block in the joint chain." % node.get_name())
                dmean = dsvc.getMean()
                if dmean > 0:
                    sn.downServiceRates[ist, r] = 1.0 / dmean

    def _refresh_local_vars(self) -> None:
        """`sn.nvars`: the width of each node's LOCAL-VARIABLE block.

        A state row is [buffer | server | vars], sliced FROM THE RIGHT, so every
        reader of a row has to know how many trailing columns are local
        variables. `sn.nvars` is that width, and this port carried the field but
        never filled it -- it stayed None, every reader took V = 0, and a row
        with local variables was decoded one block off.

        WHAT IT COST. A POLLING station's controller occupies three columns
        (position, switchover phase, k-limited counter; `polling_info(...).width`
        is the single definition). `State.toMarginal` therefore read the
        controller's own values as job counts: on polling_klimited the row
        `[0 0 0 0 2 1 0]` -- an EMPTY station whose controller is at position 2
        -- decoded as one class-1 job present, so the JMT writer emitted
        `<classPopulation population="1">` and preloaded a job that is not
        there. The MATLAB row reported QLen 1.3909 and the [M2P] row 1.5025 for
        the same model at the same seed. It only bites once a state EXISTS,
        which is why it was invisible until a solver ran before JMT and left one
        on the model.

        Mirrors `@MNetwork/refreshLocalVars.m`, whose polling arm is
        `nvars(ind, 2R+1) = State.pollingInfo(sn, ind).width`. Only the polling
        width is derived here; the MAP-restart and cache columns of the
        reference belong to their own refreshes and are added when their
        readers need them.
        """
        sn = self._sn
        R = int(sn.nclasses)
        nnodes = int(sn.nnodes)
        # The reference's shape: 3R + 1 columns, the last one shared by the
        # controllers that cannot coexist (polling, BAS blocking, breakdown).
        nvars = np.zeros((nnodes, 3 * R + 1), dtype=int)
        for ind in range(nnodes):
            if sn.isstation is None or ind >= len(sn.isstation) or not sn.isstation[ind]:
                continue
            ist = int(sn.nodeToStation[ind]) if sn.nodeToStation is not None else -1
            if ist < 0 or sn.sched is None or ist >= len(sn.sched):
                continue
            sched = sn.sched[ist]
            sched = int(sched.value) if hasattr(sched, 'value') else int(sched)
            if sched != int(SchedStrategy.POLLING.value):
                continue
            from ..api.state.polling import polling_info
            try:
                pinfo = polling_info(sn, ind)
            except Exception:
                continue  # a controller this port cannot describe declares no width
            width = pinfo.get('width') if isinstance(pinfo, dict) else getattr(pinfo, 'width', 0)
            nvars[ind, 2 * R] = int(width or 0)
        sn.nvars = nvars
        self._widen_polling_state(nvars)

    def _widen_polling_state(self, nvars) -> None:
        """Give a polling station's DEFAULT state row its controller columns.

        `sn.state` is a full local row, `[buffer | server | vars]`, and the
        default this port builds for a polling station stops after the server
        block. That was self-consistent only while `sn.nvars` was None and every
        reader assumed no local variables; now that the width is declared, a row
        that omits the block decodes to nothing at all (`toMarginal` slices the
        three controller columns off a four-column row and finds no server
        block left).

        The controller's initial value is `polling_init`, the same row
        `polling_blocks` enumerates, so the state stays a member of the space the
        generator builds. A row that already carries the block -- one that came
        from a solver, or across the JSON bridge from the reference -- is left
        alone.
        """
        sn = self._sn
        if sn.state is None or sn.nodeToStateful is None:
            return
        R = int(sn.nclasses)
        from ..api.state.polling import polling_init
        for ind in range(int(sn.nnodes)):
            if int(nvars[ind, 2 * R]) <= 0:
                continue
            isf = int(np.asarray(sn.nodeToStateful).flatten()[ind])
            if isf < 0 or isf >= len(sn.state) or sn.state[isf] is None:
                continue
            row = np.ravel(np.asarray(sn.state[isf], dtype=float))
            ist = int(sn.nodeToStation[ind])
            srvw = int(np.sum(np.asarray(sn.phasessz)[ist])) if sn.phasessz is not None else R
            want = R + srvw + int(nvars[ind, 2 * R])
            if len(row) >= want:
                continue  # already carries the block
            nbuf = row[:R] if len(row) >= R else np.zeros(R)
            try:
                ctrl = np.ravel(np.asarray(polling_init(sn, ind, nbuf, -1), dtype=float))
            except Exception:
                continue  # a controller this port cannot initialize keeps the short row
            if ctrl.size == 0:
                continue
            sn.state[isf] = np.concatenate([row, ctrl])

    def _refresh_varsparam(self) -> None:
        """Initialize sn.varsparam for cache item state tracking. Mirrors JAR
        Network.java:4745-4747 (unconditional init). -1 indicates no specific
        item selected."""
        self._sn.varsparam = -np.ones((self._sn.nnodes, 1), dtype=int)

        # Marked (MMAP) source arrivals: mark index (1-based) of class r at
        # source station i; -1 = not a marked class (mirrors MATLAB sn.markidx)
        self._sn.markidx = -np.ones((self._sn.nstations, self._sn.nclasses), dtype=int)
        any_marked = False
        for i, station in enumerate(self._stations):
            marked = getattr(station, '_marked_classes', None)
            if marked:
                any_marked = True
                for k, cls in enumerate(marked):
                    j = self._classes.index(cls)
                    self._sn.markidx[i, j] = k + 1

        # Marked (MMAP) source classes share the carrier's modulating chain, so a
        # non-carrier mark contributes a single always-zero state column rather
        # than its own phase block (refreshProcessRepresentations.m). Applied
        # here because markidx is only known after _refresh_rates has run; without
        # it the carrier's phase never enters the state space and every marked
        # class reports the UNWEIGHTED sum of its per-phase arrival rates.
        if any_marked and self._sn.phasessz is not None:
            self._sn.phasessz = np.asarray(self._sn.phasessz).astype(int)
            self._sn.phasessz[self._sn.markidx > 1] = 1
            self._sn.phaseshift = np.hstack([
                np.zeros((self._sn.nstations, 1), dtype=int),
                np.cumsum(self._sn.phasessz, axis=1)
            ]).astype(int)

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
            sn.regionsz = np.ones((F, K), dtype=float)
            # regionlincon[f] = [A, b] linear constraint pair on the same row for region f
            sn.regionlincon = [None] * F
            sn.regionmaxmem = []
            # region membership recorded explicitly (region[f]==-1 is ambiguous); see _kb/04-networkstruct.md Python native network.py section.
            sn.regionmembers = []

            for f, fcr in enumerate(self._regions):
                region_matrix = -1 * np.ones((M, K + 1), dtype=float)
                region_mem_matrix = -1 * np.ones((M, 1), dtype=float)
                region_member_mask = np.zeros(M, dtype=bool)

                for node in fcr.nodes:
                    for i, station in enumerate(self._stations):
                        if station is node:
                            region_member_mask[i] = True
                            # per-class memory budget folds to a job cap floor(maxMem_r/classSize_r); JMT classMemoryConstraint semantics.
                            for r, job_class in enumerate(self._classes):
                                cap_r = fcr.get_class_max_jobs(job_class)
                                mem_r = fcr.get_class_max_memory(job_class)
                                sz_r = fcr.get_class_size(job_class)
                                if mem_r != -1 and sz_r > 0:
                                    memjobs = int(mem_r // sz_r)
                                    cap_r = memjobs if cap_r == -1 else min(cap_r, memjobs)
                                region_matrix[i, r] = cap_r
                            region_matrix[i, K] = fcr.global_max_jobs
                            # Replicate the region-global memory budget on this member row
                            region_mem_matrix[i, 0] = fcr.global_max_memory
                            break

                for r, job_class in enumerate(self._classes):
                    drop_rule = fcr.get_drop_rule(job_class)
                    sn.regionrule[f, r] = float(drop_rule)
                    sn.regionweight[f, r] = fcr.get_class_weight(job_class)
                    sn.regionsz[f, r] = fcr.get_class_size(job_class)

                sn.region.append(region_matrix)
                sn.regionmaxmem.append(region_mem_matrix)
                sn.regionmembers.append(region_member_mask)

                # Serialize the linear constraint pair (A,b) on the same row if set
                if fcr.has_linear_constraints():
                    linConA, linConB = fcr.get_linear_constraints()
                    sn.regionlincon[f] = [linConA, linConB]
        else:
            sn.nregions = 0
            sn.region = []
            sn.regionrule = np.array([])
            sn.regionweight = np.array([])
            sn.regionsz = np.array([])
            sn.regionmaxmem = []
            sn.regionmembers = []

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

        # Node types. A QUEUE WITH INFINITE SERVERS IS A DELAY STATION, as
        # MATLAB's getNodeTypes and the JAR's Network.getNodeTypes both report
        # it: the analyzers partition the stations on this field, and a station
        # that serves every job on arrival belongs with the delays whichever
        # class built it.
        from .base import NodeType as _NT
        self._sn.nodetype = [
            _NT.DELAY if (node.node_type == _NT.QUEUE
                          and np.isinf(getattr(node, 'number_of_servers', 1)))
            else node.node_type
            for node in self._nodes]

        # Station classification
        self._sn.isstation = np.array([node.is_station() for node in self._nodes], dtype=bool)

        # Stateful classification (nodes that maintain state)
        # Stateful: SOURCE, DELAY, QUEUE, CACHE, JOIN, ROUTER, PLACE, TRANSITION
        from .base import NodeType
        stateful_types = {NodeType.SOURCE, NodeType.DELAY, NodeType.QUEUE,
                         NodeType.CACHE, NodeType.JOIN, NodeType.ROUTER,
                         NodeType.PLACE, NodeType.TRANSITION}
        # FJ tag-augmented copies treat Fork nodes as stateful (they hold the
        # parent job for one vanishing state before the fork firing)
        if getattr(self, 'fork_stateful', False):
            stateful_types = stateful_types | {NodeType.FORK}
        self._sn.isstateful = np.array([node.node_type in stateful_types for node in self._nodes], dtype=bool)
        self._sn.nstateful = int(np.sum(self._sn.isstateful))
        # Station-indexed mask of queues carrying setup/delay-off times, as in
        # MATLAB refreshStruct.m.
        from .nodes import Queue
        self._sn.hassetup = np.array(
            [isinstance(station, Queue) and bool(getattr(station, '_setup_time', None))
             for station in self._stations], dtype=bool)
        # FJ tag-augmented copy marker (see ModelAdapter.fjtag): Join/Fork
        # carry per-class count-vector states
        self._sn.isfjaugmented = bool(getattr(self, 'is_fj_augmented', False))

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
                is_open_class = (jobclass.jobclass_type == JobClassType.OPEN)
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

            # open classes must be detected explicitly: JobClass.getNumberOfJobs() returns 0 for them, which would otherwise collapse classcap/chaincap to 0.
            if jobclass.jobclass_type == JobClassType.OPEN:
                njobs[j] = np.inf
            elif hasattr(jobclass, '_njobs'):
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

        # Self-looping class flag (mirrors MATLAB refreshStruct): true when the
        # class perpetually cycles at its reference station (SelfLoopingClass).
        isslc = np.zeros(nclasses, dtype=bool)
        for j, jobclass in enumerate(self._classes):
            if jobclass.__class__.__name__ == 'SelfLoopingClass':
                isslc[j] = True
        self._sn.isslc = isslc

        # signaltarget(c) = 0-based index of the (positive) class whose jobs signal
        # class c removes; -1 for non-signals or signals with no target.
        signaltarget = np.full(nclasses, -1, dtype=int)
        for j, jobclass in enumerate(self._classes):
            if issignal[j] and hasattr(jobclass, 'getTargetJobClass'):
                tgt = jobclass.getTargetJobClass()
                if tgt is not None:
                    for jj, cc in enumerate(self._classes):
                        if cc is tgt:
                            signaltarget[j] = jj
                            break
        self._sn.signaltarget = signaltarget

        # Signal removal properties
        signalremdist = [None] * nclasses
        signalrempolicy = [None] * nclasses
        iscatastrophe = np.zeros(nclasses, dtype=bool)
        for j, jobclass in enumerate(self._classes):
            if issignal[j]:
                if hasattr(jobclass, '_removal_distribution'):
                    signalremdist[j] = jobclass._removal_distribution
                if hasattr(jobclass, '_removal_policy'):
                    signalrempolicy[j] = jobclass._removal_policy
                if hasattr(jobclass, '_signal_type'):
                    from .classes import SignalType
                    if jobclass._signal_type == SignalType.CATASTROPHE:
                        iscatastrophe[j] = True
        self._sn.signalremdist = signalremdist
        self._sn.signalrempolicy = signalrempolicy
        self._sn.iscatastrophe = iscatastrophe

        # Immediate feedback matrix (station x class)
        nstations = len(self._stations)
        immfeed = np.zeros((nstations, nclasses), dtype=bool)
        from .nodes import Queue
        for ist, station in enumerate(self._stations):
            node = station
            for r, jobclass in enumerate(self._classes):
                station_has = False
                if isinstance(node, Queue) and hasattr(node, '_immediate_feedback_classes'):
                    if node._immediate_feedback_all:
                        station_has = True
                    else:
                        station_has = r in node._immediate_feedback_classes
                class_has = jobclass._immediate_feedback
                immfeed[ist, r] = station_has or class_has
        self._sn.immfeed = immfeed

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

        # True where the representation is Markovian, so that mu, phi and pie
        # carry their probabilistic reading. Consumers of those fields (CTMC
        # state space, SSA, fluid ODEs) treat them as rates and probabilities,
        # which fails for a matrix-exponential process: its per-phase values are
        # signed. Filled per station-class by _extract_process_params.
        self._sn.isph = np.ones((nstations, nclasses), dtype=bool)

        # Process type IDs: (M x K) array of ProcessType values
        self._sn.procid = np.empty((nstations, nclasses), dtype=object)
        self._sn.procid.fill(ProcessType.EXP)  # Default to exponential

        # phases/phasessz mirror MATLAB refreshProcessPhases; sn.phases was previously never populated here.
        self._sn.phases = np.zeros((nstations, nclasses))

        # per (station,class) Laplace-Stieltjes transforms for the G/M/1 MVA sigma-root; mirrors MATLAB refreshLST/JAR sn.lst.
        self._sn.lst = [[None for _ in range(nclasses)] for _ in range(nstations)]

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
                    # phase count from the actual distribution, not re-derived from SCV; see _kb/04-networkstruct.md Python native network.py section.
                    self._sn.phases[i, j] = self._dist_num_phases(dist)

                    # Immediate distribution (mean=0) uses GlobalConstants.Immediate rather than inf, matching MATLAB.
                    if hasattr(dist, 'isImmediate') and dist.isImmediate():
                        self._sn.rates[i, j] = GlobalConstants.Immediate
                    elif hasattr(dist, 'isDisabled') and dist.isDisabled():
                        # Disabled distribution: rate=NaN (not 1/inf=0), matching MATLAB's NaN handling.
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

                    # marked (MMAP) arrival: rate/SCV come from the mark's marginal MAP, mirrors MATLAB refreshRates.
                    from ..distributions.markovian import MarkedMAP as _MarkedMAP
                    _marked = getattr(station, '_marked_classes', None)
                    if isinstance(dist, _MarkedMAP) and _marked and jobclass in _marked:
                        _marginal = dist.to_marginal_map(_marked.index(jobclass) + 1)
                        self._sn.rates[i, j] = _marginal.getRate()
                        self._sn.scv[i, j] = _marginal.getSCV()

                    # a class already flagged DISABLED keeps that procid; MATLAB refreshProcessTypes likewise never resolves a type for it.
                    if not (hasattr(dist, 'isDisabled') and dist.isDisabled()):
                        self._extract_process_params(i, j, dist)
                        # The extractor returns from many branches, so the
                        # Markovian test is applied here, where every branch has
                        # already written sn.proc for this station-class.
                        from ..api.sn.utils import sn_is_phasetype
                        self._sn.isph[i, j] = sn_is_phasetype(self._sn.proc[i][j])

                    # Laplace-Stieltjes transform closure for the G/M/1 sigma-root
                    # (non-PH arrivals). Bind dist by default arg to capture it.
                    if hasattr(dist, 'evalLST'):
                        self._sn.lst[i][j] = (lambda s, _d=dist: _d.evalLST(s))
                else:
                    # no distribution: Source has rate=0 (no arrival), Join/Cache have no/instant service work.
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

        # phasessz/phaseshift: disabled class still takes one slot; mirrors MATLAB refreshProcessPhases.
        self._sn.phasessz = np.maximum(self._sn.phases, 1).astype(int)
        # A Join keeps its true phase count, including the zero of a class it
        # never serves; mirrors refreshProcessRepresentations.m.
        if self._sn.nodetype is not None:
            for ind in range(self._sn.nnodes):
                if self._sn.nodetype[ind] != NodeType.JOIN:
                    continue
                ist = int(self._sn.nodeToStation[ind])
                if 0 <= ist < nstations:
                    self._sn.phasessz[ist, :] = self._sn.phases[ist, :].astype(int)
        self._sn.phaseshift = np.hstack([
            np.zeros((nstations, 1), dtype=int),
            np.cumsum(self._sn.phasessz, axis=1)
        ]).astype(int)

    @staticmethod
    def _dist_num_phases(dist) -> int:
        """Number of phases of a distribution, 0 when it is disabled.

        Mirrors MATLAB refreshProcessPhases, which records length(mu_i{r}) and
        0 for a disabled process.
        """
        try:
            if hasattr(dist, 'isDisabled') and dist.isDisabled():
                return 0
        except Exception:
            pass
        getn = getattr(dist, 'getNumberOfPhases', None)
        if getn is not None:
            try:
                return int(getn())
            except Exception:
                pass
        return 1

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
        from ..distributions.markovian import MAP, MMPP2, PH, APH, Coxian, Cox2, BMAP, MarkedMAP, ME, RAP
        from ..distributions.continuous import (
            Exp, Erlang, HyperExp, Gamma, Uniform, Det, Pareto, Weibull, Lognormal, Immediate,
            NHPP, MAPt, PHt
        )
        from ..distributions.discrete import Geometric, Bernoulli, Binomial, Poisson

        # Check for BMAP first (batch Markovian; extends MarkedMAP, not MAP)
        if isinstance(dist, BMAP):
            self._sn.procid[station_idx, class_idx] = ProcessType.BMAP
            # Layout mirrors the JAR MatrixCell: [D0, D1, D_batch1, ..., D_batchK]
            batch_mats = [np.atleast_2d(np.asarray(Dk, dtype=float)) for Dk in dist._process[1:]]
            D0 = np.atleast_2d(np.asarray(dist._process[0], dtype=float))
            D1 = np.sum(batch_mats, axis=0)
            self._sn.proc[station_idx][class_idx] = [D0, D1] + batch_mats

        # marked MAP (MMAP): shared modulating chain, M3A layout [D0, D1_agg, D11,...,D1K].
        elif isinstance(dist, MarkedMAP):
            self._sn.procid[station_idx, class_idx] = ProcessType.MMAP
            self._sn.proc[station_idx][class_idx] = dist.to_m3a()

        # Discrete-time (D0,D1): the pair already has MAP shape, so it reaches
        # sn.proc verbatim. Without this branch a DMAP fell through to the
        # Replayer catch-all and reached sn.proc as None, while sn.rates still
        # looked right. It must precede MAP even though DMAP does not subclass
        # it, so the reading stays symmetric with the MATLAB and JAR arms.
        elif type(dist).__name__ == 'DMAP':
            self._sn.procid[station_idx, class_idx] = ProcessType.DMAP
            self._sn.proc[station_idx][class_idx] = [
                np.atleast_2d(np.asarray(dist.getD0(), dtype=float)),
                np.atleast_2d(np.asarray(dist.getD1(), dtype=float))]

        # Check for MMPP2 first (more specific than MAP)
        elif isinstance(dist, MMPP2):
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

        # RAP/ME carry a (D0,D1) pair consumed directly by the map_* algorithms and RAP/RAP/1 QBD.
        elif isinstance(dist, RAP):
            self._sn.procid[station_idx, class_idx] = ProcessType.RAP
            self._sn.proc[station_idx][class_idx] = [dist.getD0(), dist.getD1()]

        elif isinstance(dist, ME):
            self._sn.procid[station_idx, class_idx] = ProcessType.ME
            self._sn.proc[station_idx][class_idx] = [dist.getD0(), dist.getD1()]

        # Coxian/Cox2 branch ordering: see _kb/04-networkstruct.md Python native network.py section (class-hierarchy trap).
        # Stored as (D0,D1) like every other Markovian family, not as (alpha,T):
        # the two are indistinguishable by shape, so a reader that assumed the
        # pair was a MAP read alpha as D0.
        elif isinstance(dist, Cox2):
            self._sn.procid[station_idx, class_idx] = ProcessType.COX2
            alpha = dist.getInitProb() if hasattr(dist, 'getInitProb') else None
            T = dist.getD0() if hasattr(dist, 'getD0') else None
            if alpha is not None and T is not None:
                from ..api.sn.proc_form import map_from_ph
                self._sn.proc[station_idx][class_idx] = list(map_from_ph(alpha, T))

        elif isinstance(dist, Coxian):
            self._sn.procid[station_idx, class_idx] = ProcessType.COXIAN
            alpha = dist.getInitProb() if hasattr(dist, 'getInitProb') else None
            T = dist.getD0() if hasattr(dist, 'getD0') else None
            if alpha is not None and T is not None:
                from ..api.sn.proc_form import map_from_ph
                self._sn.proc[station_idx][class_idx] = list(map_from_ph(alpha, T))

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
                from ..api.sn.proc_form import map_from_ph
                self._sn.proc[station_idx][class_idx] = list(map_from_ph(alpha, T))

        # Check for HyperExp
        elif isinstance(dist, HyperExp):
            self._sn.procid[station_idx, class_idx] = ProcessType.HYPEREXP
            # Store parameters for HyperExp (uses _probs and _rates arrays)
            if hasattr(dist, '_probs') and hasattr(dist, '_rates'):
                from ..api.sn.proc_form import map_from_hyperexp
                self._sn.proc[station_idx][class_idx] = list(
                    map_from_hyperexp(dist._probs, dist._rates))

        # Check for Erlang
        elif isinstance(dist, Erlang):
            self._sn.procid[station_idx, class_idx] = ProcessType.ERLANG
            # Erlang uses _phases and _phase_rate attributes
            if hasattr(dist, '_phases') and hasattr(dist, '_phase_rate'):
                from ..api.sn.proc_form import map_from_erlang
                self._sn.proc[station_idx][class_idx] = list(
                    map_from_erlang(dist._phases, dist._phase_rate))

        # Check for Exponential
        elif isinstance(dist, Exp):
            self._sn.procid[station_idx, class_idx] = ProcessType.EXP
            # For Exp, proc can remain None or store rate
            if hasattr(dist, '_rate'):
                from ..api.sn.proc_form import map_from_exp
                self._sn.proc[station_idx][class_idx] = list(map_from_exp(dist._rate))

        # non-Markovian distributions: raw params stored in sn.proc (MATLAB layout) so sn_nonmarkov_toph can refit the PDF.
        elif isinstance(dist, Gamma):
            self._sn.procid[station_idx, class_idx] = ProcessType.GAMMA
            self._sn.proc[station_idx][class_idx] = [dist._shape, dist._scale]
        elif isinstance(dist, Uniform):
            self._sn.procid[station_idx, class_idx] = ProcessType.UNIFORM
            self._sn.proc[station_idx][class_idx] = [dist._min, dist._max]
        elif isinstance(dist, Immediate):
            # Check Immediate before Det since Immediate extends Det
            self._sn.procid[station_idx, class_idx] = ProcessType.IMMEDIATE
        elif isinstance(dist, Det):
            self._sn.procid[station_idx, class_idx] = ProcessType.DET
            self._sn.proc[station_idx][class_idx] = [dist._value]
        elif isinstance(dist, Pareto):
            self._sn.procid[station_idx, class_idx] = ProcessType.PARETO
            self._sn.proc[station_idx][class_idx] = [dist._alpha, dist._scale]
        elif isinstance(dist, Weibull):
            self._sn.procid[station_idx, class_idx] = ProcessType.WEIBULL
            self._sn.proc[station_idx][class_idx] = [dist._shape, dist._scale]
        elif isinstance(dist, Lognormal):
            self._sn.procid[station_idx, class_idx] = ProcessType.LOGNORMAL
            self._sn.proc[station_idx][class_idx] = [dist._mu, dist._sigma]
        elif isinstance(dist, (Bernoulli, Binomial, Poisson)):
            # counting-distribution interval: {mean,SCV} pair; LDES resolves the zero atom to an immediate interval.
            self._sn.procid[station_idx, class_idx] = (
                ProcessType.BERNOULLI if isinstance(dist, Bernoulli)
                else ProcessType.BINOMIAL if isinstance(dist, Binomial)
                else ProcessType.POISSON)
            self._sn.proc[station_idx][class_idx] = [
                dist.getMean(), dist.getVar() / (dist.getMean() ** 2)]
        elif isinstance(dist, Geometric):
            # Geometric-lattice interarrival/service: {mean,SCV} pair, no MAP representation, matching MATLAB/JAR.
            self._sn.procid[station_idx, class_idx] = ProcessType.GEOMETRIC
            self._sn.proc[station_idx][class_idx] = [dist.getMean(), dist.getVar() / (dist.getMean() ** 2)]
        elif isinstance(dist, NHPP):
            self._sn.procid[station_idx, class_idx] = ProcessType.NHPP
            # rate schedule stored as {breakpoints, rates, cyclic}, matching MATLAB cell / JAR MatrixCell.
            self._sn.proc[station_idx][class_idx] = [
                np.atleast_2d(np.asarray(dist.breakpoints, dtype=float)),
                np.atleast_2d(np.asarray(dist.rates, dtype=float)),
                bool(dist.cyclic),
            ]
        elif isinstance(dist, MAPt):
            self._sn.procid[station_idx, class_idx] = ProcessType.MAPT
            # flat layout [breakpoints, D0_1..D0_n, D1_1..D1_n, cyclic], so the
            # JAR MatrixCell (which cannot nest) carries the same slot; n follows
            # from (len-2)//2.
            self._sn.proc[station_idx][class_idx] = (
                [np.atleast_2d(np.asarray(dist.breakpoints, dtype=float))]
                + [np.asarray(M, dtype=float) for M in dist.D0]
                + [np.asarray(M, dtype=float) for M in dist.D1]
                + [bool(dist.cyclic)])
        elif isinstance(dist, PHt):
            self._sn.procid[station_idx, class_idx] = ProcessType.PHT
            # flat layout [breakpoints, alpha_1..alpha_n, S_1..S_n, cyclic], each
            # alpha as a 1-by-h row so every element of the cell is a matrix.
            self._sn.proc[station_idx][class_idx] = (
                [np.atleast_2d(np.asarray(dist.breakpoints, dtype=float))]
                + [np.atleast_2d(np.asarray(a, dtype=float)) for a in dist.alpha]
                + [np.asarray(M, dtype=float) for M in dist.S]
                + [bool(dist.cyclic)])

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
            if proc_type is None:
                # unrecognized ProcessType.fromText mirrors MATLAB's error rather than silently defaulting to EXP.
                raise ValueError('Unrecognized process type: %s' % class_name)
            self._sn.procid[station_idx, class_idx] = proc_type

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
        elif isinstance(sched, int):
            # a raw native enum value, e.g. from LayeredNetworkStruct.sched
            return SchedStrategy(sched)
        else:
            name = str(sched).split('.')[-1]

        # Case-insensitive: callers may pass a lowercase/mixed-case string
        # (e.g. Queue(model, name, "fcfs")) alongside a canonical enum/name;
        # every SchedStrategy member name is uppercase, so this only widens
        # what's accepted -- it never changes the outcome for already-correct
        # casing.
        name = name.upper()

        # OI is a PAS specialization with an empty swap graph; normalized to PAS so all dispatch sites share one code path.
        if name == 'OI':
            name = 'PAS'

        # FCFSPRIO normalized to HOL (same enum value in MATLAB) so dispatch sites testing sched==HOL still match.
        if name == 'FCFSPRIO':
            name = 'HOL'

        # Map to native SchedStrategy by name
        try:
            return SchedStrategy[name]
        except KeyError:
            raise ValueError("Unknown scheduling strategy %r; it cannot be normalized to a "
                             "native SchedStrategy." % (sched,))

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
        # Cox2/Coxian/PH must be tested before the ME/RAP class-name fallback.
        from ..distributions.markovian import APH, PH, MAP, MMPP2, Cox2

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
        elif isinstance(dist, Cox2):
            return ProcessType.COX2
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
            if proc_type is None:
                # Mirrors MATLAB ProcessType.fromText: an unrecognized process
                # type is an error, not an exponential. See _extract_process_params.
                raise ValueError('Unrecognized process type: %s' % class_name)
            return proc_type

    def _refresh_load_dependence(self) -> None:
        """Extract load-dependent scaling from stations.

        Mirrors MATLAB getLimitedLoadDependence.m: the lattice width is the
        LONGEST user-supplied lldScaling vector, shorter rows stay at 1, and
        zeros are replaced by ones. It does NOT depend on the closed population:
        sizing by sum(njobs) instead dropped the handle entirely for all-OPEN
        models (total_pop == 0), silently returning the load-independent answer.
        """
        nstations = len(self._stations)

        _as_1d = lld_scaling_as_1d

        # Lattice width = longest configured scaling vector (MATLAB maxsize)
        maxsize = 0
        for station in self._stations:
            ld = station.get_load_dependence() if hasattr(station, 'get_load_dependence') else None
            if ld is not None and len(ld) > 0:
                maxsize = max(maxsize, len(_as_1d(ld)))

        if maxsize == 0:
            # no load dependence: None is the python spelling of MATLAB's ones(M,0) (size(.,2)==0 sentinel).
            self._sn.lldscaling = None
            return

        lldscaling = np.ones((nstations, maxsize))
        for i, station in enumerate(self._stations):
            ld = station.get_load_dependence() if hasattr(station, 'get_load_dependence') else None
            if ld is not None and len(ld) > 0:
                ld_array = _as_1d(ld)
                lldscaling[i, :len(ld_array)] = ld_array
                lldscaling[i, lldscaling[i, :] == 0] = 1.0

        self._sn.lldscaling = lldscaling

    def _refresh_class_dependence(self) -> None:
        """Extract class-dependence handles beta_{i,r}(n) from stations.

        Stations that declare no class dependence are deliberately left as None
        rather than filled with a constant 1: the entry is a class-dependent
        RATE (Sauer 1983, eq. (40)), for which a constant 1 would assert "every
        class completes at rate 1" -- not load independence (an LI station is
        beta_{i,r}(n) = mu_i * n_r/|n|, whose n_r/|n| factor is what regenerates
        the multinomial). Consumers treat a None entry as "no class dependence":
        pfqn_cdfun skips it and pfqn_conv keeps the station on the
        load-independent recurrence. cdscaling is set to None entirely when no
        station declares one, so a falsy cdscaling still means "no class
        dependence anywhere".
        """
        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Check if any station has class dependence
        has_cd = False
        cdscaling_list = [None] * nstations
        # Declared peak rate scaling per class (Util = T*S/peak); NaN where a
        # station has no class dependence (broadcast a scalar to all classes).
        cdscalingpeak = np.full((nstations, nclasses), np.nan)

        for i, station in enumerate(self._stations):
            cd = station.get_class_dependence() if hasattr(station, 'get_class_dependence') else None
            if cd is not None:
                has_cd = True
                cdscaling_list[i] = cd
                pk = station.get_class_dependence_peak() if hasattr(station, 'get_class_dependence_peak') else None
                if pk is None:
                    raise ValueError(
                        "Class-dependent station %d has no declared peak rate; "
                        "use set_class_dependence(beta, peak_rate_per_class)." % i)
                pk = np.atleast_1d(np.asarray(pk, dtype=float)).ravel()
                if pk.size == 1:
                    cdscalingpeak[i, :] = pk[0]
                elif pk.size == nclasses:
                    cdscalingpeak[i, :] = pk
                else:
                    raise ValueError(
                        "peak_rate_per_class at station %d must be scalar or of "
                        "length nclasses=%d." % (i, nclasses))

        if has_cd:
            self._sn.cdscaling = cdscaling_list
            self._sn.cdscalingpeak = cdscalingpeak
        else:
            self._sn.cdscaling = None
            self._sn.cdscalingpeak = None

        # Joint-dependence handles eta_i(n) (non-product-form). Same extraction
        # as class dependence; kept in a separate field/peak to preserve the
        # product-form vs joint distinction (see set_joint_dependence).
        has_jd = False
        jdscaling_list = [None] * nstations
        jdscalingpeak = np.full((nstations, nclasses), np.nan)
        for i, station in enumerate(self._stations):
            jd = station.get_joint_dependence() if hasattr(station, 'get_joint_dependence') else None
            if jd is not None:
                has_jd = True
                jdscaling_list[i] = jd
                pk = station.get_joint_dependence_peak() if hasattr(station, 'get_joint_dependence_peak') else None
                if pk is None:
                    raise ValueError(
                        "Joint-dependent station %d has no declared peak rate; "
                        "use set_joint_dependence(eta, peak_rate_per_class)." % i)
                pk = np.atleast_1d(np.asarray(pk, dtype=float)).ravel()
                if pk.size == 1:
                    jdscalingpeak[i, :] = pk[0]
                elif pk.size == nclasses:
                    jdscalingpeak[i, :] = pk
                else:
                    raise ValueError(
                        "peak_rate_per_class at station %d must be scalar or of "
                        "length nclasses=%d." % (i, nclasses))

        if has_jd:
            self._sn.jdscaling = jdscaling_list
            self._sn.jdscalingpeak = jdscalingpeak
        else:
            self._sn.jdscaling = None
            self._sn.jdscalingpeak = None

        # Network-level global dependence phi(n) over the FULL population matrix
        # (the Whittle primitive); see set_global_dependence.
        self._sn.gdscaling = self._gd_scaling
        self._sn.gdscalingpeak = self.get_global_dependence_peak()
        self._sn.gdscalingcutoff = self.get_global_dependence_cutoff()

    def set_global_dependence(self, phi, peak, cutoff=10) -> None:
        """Declare a globally state-dependent service-rate scaling phi(n).

        n is the FULL (nstations, nclasses) population matrix, not the population
        local to one station. This is the Whittle-network primitive: when phi
        satisfies phi_s(n) phi_t(n-e_s) = phi_t(n) phi_s(n-e_t) the chain is
        reversible with pi(n) ~ Phi(n) prod rho_s**n_s and is insensitive. It also
        expresses bandwidth sharing, where one route holds several links at once
        and no per-station scaling can reproduce the coupling.

        phi returns a scalar (broadcast), an (nstations,) column (per station) or
        an (nstations, nclasses) matrix. The effective rate of class r at station i
        is its base rate times phi[i, r], composing multiplicatively with any
        load-, class- or joint-dependence. Only SolverCTMC supports it.

        Args:
            phi: callable over the full population matrix
            peak: REQUIRED peak scaling (scalar, (nstations,) or (nstations, nclasses)),
                  normalizing utilization as Util = T*S/peak
            cutoff: per-slot OPEN-class truncation used when phi is materialized
                  onto the JSON wire (closed classes are tabulated up to their own
                  population). It plays no part in solving, and exists because a
                  handle cannot cross a language boundary: the writer needs to know
                  how far the lattice extends. Set it to the cutoff the model is
                  solved at.
        """
        if not callable(phi):
            raise ValueError("Global dependence must be specified through a callable.")
        if peak is None:
            raise ValueError(
                "Global dependence requires an explicit peak rate: "
                "set_global_dependence(phi, peak).")
        peak_arr = np.atleast_1d(np.asarray(peak, dtype=float))
        if np.any(peak_arr <= 0):
            raise ValueError("peak must be positive.")
        M = len(self._stations)
        K = len(self._classes)
        # Probe now so a wrong output shape is refused at declaration time rather
        # than midway through state-space generation.
        for probe in (np.zeros((M, K)), np.ones((M, K))):
            v = np.asarray(phi(probe), dtype=float)
            if not np.all(np.isfinite(v)) or np.any(v < 0):
                raise ValueError(
                    "The global dependence handle must return finite nonnegative scalings.")
            if v.size != 1 and v.shape != (M,) and v.shape != (M, 1) and v.shape != (M, K):
                raise ValueError(
                    "The global dependence handle must return a scalar, an (%d,) column "
                    "or an (%d, %d) matrix." % (M, M, K))
        if not isinstance(cutoff, (int, np.integer)) or int(cutoff) < 1:
            raise ValueError("cutoff must be a positive integer.")
        self._gd_scaling = phi
        self._gd_scaling_peak = peak_arr
        self._gd_scaling_cutoff = int(cutoff)
        self.reset_struct()

    def get_global_dependence(self):
        """Network-level global dependence handle, or None if the model declares none."""
        return self._gd_scaling

    def get_global_dependence_cutoff(self):
        """Per-slot open-class wire truncation of the global dependence, or None if none is declared."""
        if self._gd_scaling is None:
            return None
        return int(getattr(self, '_gd_scaling_cutoff', 10) or 10)

    def get_global_dependence_peak(self):
        """(nstations, nclasses) peak of the global dependence, or None if none is declared."""
        if self._gd_scaling is None:
            return None
        M = len(self._stations)
        K = len(self._classes)
        pk = np.asarray(self._gd_scaling_peak, dtype=float)
        out = np.ones((M, K))
        if pk.size == 1:
            out[:, :] = float(pk.ravel()[0])
        elif pk.shape == (M,) or pk.shape == (M, 1):
            out[:, :] = pk.reshape(M, 1)
        elif pk.shape == (M, K):
            out[:, :] = pk
        else:
            raise ValueError(
                "peak must be a scalar, an (%d,) column or an (%d, %d) matrix." % (M, M, K))
        return out

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

            # LPS concurrency limit in schedparam column 0, mirrors MATLAB Queue.setLimit; state engine still caps at nservers.
            if self._sn.sched[i] == SchedStrategy.LPS.value and \
                    getattr(station, '_lps_limit', None) is not None:
                self._sn.schedparam[i, 0] = station._lps_limit

    def _quorum_join_classes(self) -> set:
        """Class indices whose Join fires on a STRICT quorum, i.e. on fewer siblings
        than are forked. Such a class is not population-conserving at the sibling level.

        Siblings are counted at the FORK, as the engines count them: its out-degree
        times tasksPerLink, with the join in-degree as the fallback.
        """
        from .nodes import Join as _JoinNode, Fork as _ForkNode
        from .base import JoinStrategy as _JS

        out = set()
        if not self._nodes:
            return out
        conn = self.get_connection_matrix()
        for node_idx, node in enumerate(self._nodes):
            if not isinstance(node, _JoinNode):
                continue
            nsib = int(np.count_nonzero(conn[:, node_idx]))
            fork = node.get_fork()
            if isinstance(fork, _ForkNode):
                fork_idx = self._nodes.index(fork)
                tpl = fork.get_tasks_per_link()
                w = 1
                if tpl is not None:
                    tpl = np.atleast_1d(tpl)
                    if tpl.size > 0:
                        w = max(1, int(round(float(tpl.flat[0]))))
                nsib = int(np.count_nonzero(conn[fork_idx, :])) * w
            for r, jobclass in enumerate(self._classes):
                if node.get_strategy(jobclass) == _JS.STD:
                    continue
                kreq = node.get_required(jobclass)
                if kreq is None or kreq <= 0:
                    continue
                if nsib <= 0 or kreq < nsib:
                    out.add(int(r))
        return out

    def _fork_tasks_per_link_factor(self) -> int:
        """Product of tasksPerLink over every Fork, the factor by which a fork can
        multiply the jobs a branch station holds per circulating parent.

        The PRODUCT rather than the max, because a fork nested in another's branch
        multiplies again; for forks in series it over-bounds, and a cap that never
        binds costs nothing.
        """
        from .nodes import Fork as _ForkNode

        factor = 1
        for node in self._nodes:
            if not isinstance(node, _ForkNode):
                continue
            tpl = node.get_tasks_per_link()
            if tpl is None:
                continue
            tpl = np.atleast_1d(np.asarray(tpl, dtype=float))
            if tpl.size > 0 and np.isfinite(tpl.flat[0]):
                factor *= max(1, int(round(float(tpl.flat[0]))))
        return factor

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

        # a chain routed through a QUORUM join is not population-conserving either, and for the
        # same reason as a spawn-target chain: the join releases the parent at the k-th of n
        # siblings and the n-k stragglers stay in the branches, so the parent forks again while
        # they are still in flight. Nothing bounds that backlog, so a branch station holds no
        # more than the class population only under a STANDARD join; capping it at sum(njobs)
        # makes a simulator drop a closed job. Read off the node objects, not off sn.nodeparam:
        # _refresh_capacity runs before _refresh_local_vars rebuilds it.
        # see _kb/05-solvers-overview.md
        quorum_classes = self._quorum_join_classes()
        # a fork with tasksPerLink = w > 1 puts w tasks of the SAME parent on one link,
        # so a branch station can hold w jobs per circulating parent and the chain
        # population is no longer its bound. The multiplier is the PRODUCT over the
        # forks, because a fork nested in another's branch multiplies again; that is an
        # upper bound for forks in series, where a cap that never binds costs nothing,
        # and exact for the single-fork case. Without it a simulator drops a closed job
        # at a branch station. see _kb/04-networkstruct.md
        fork_task_factor = self._fork_tasks_per_link_factor()

        # Compute chain-based capacities
        for c in range(nchains):
            # Get classes in this chain
            if inchain is not None and c < len(inchain):
                classes_in_chain = inchain[c]
            else:
                classes_in_chain = list(range(nclasses))

            # chain capacity is the summed population, including Inf for any open class in the chain.
            chain_cap = np.sum([njobs[r] for r in classes_in_chain]) * fork_task_factor

            # a chain holding a spawn-target class is not population-conserving and must not get a finite capacity.
            for jc in self._classes:
                spawn_cls = getattr(jc, '_spawn_class', None)
                if spawn_cls is not None and spawn_cls in self._classes and \
                        self._classes.index(spawn_cls) in list(classes_in_chain):
                    chain_cap = np.inf
                    break

            if any(int(r) in quorum_classes for r in classes_in_chain):
                chain_cap = np.inf

            for r in classes_in_chain:
                for ist in range(nstations):
                    station = self._stations[ist]

                    # NaN rate means disabled, except Place nodes (no service rate); mirrors MATLAB refreshCapacity.m.
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

                        # a finite orbit bounds station population at servers+orbit capacity; see _kb/04-networkstruct.md droprule/classcap section.
                        omax_map = getattr(station, '_orbit_max_jobs', None)
                        if omax_map is not None and jobclass is not None:
                            omax = omax_map.get(jobclass, -1)
                            if omax is not None and omax >= 0:
                                nsrv = station.number_of_servers
                                if nsrv is None or not np.isfinite(nsrv):
                                    nsrv = 1
                                classcap[ist, r] = min(classcap[ist, r], int(nsrv) + int(omax))

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

        # droprule defaults: WAITQ for infinite capacity, DROP for finite (unless explicit); mirrors MATLAB refreshCapacity.m.
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
                    elif not isinstance(station._drop_rule, (dict, list)):
                        # Scalar drop rule (e.g., Station._drop_rule is a single DropStrategy)
                        station_drop_rule = station._drop_rule

                # explicit WAITQ at an open-class finite buffer is rejected (no solver honours it); see _kb/04-networkstruct.md droprule/classcap section.
                _is_user_rule_node = type(station).__name__ not in ('Place', 'Join')
                if _is_user_rule_node and station_drop_rule is not None:
                    from .base import DropStrategy as _BaseDropStrategy
                    _is_waitq = (station_drop_rule == _BaseDropStrategy.WaitingQueue
                                 or (not isinstance(station_drop_rule, _BaseDropStrategy)
                                     and int(station_drop_rule) == int(SnDropStrategy.WAITQ)))
                    _njobs_r = float(np.asarray(self._sn.njobs, dtype=float).ravel()[r]) \
                        if getattr(self._sn, 'njobs', None) is not None else np.inf
                    if _is_waitq and np.isinf(_njobs_r):
                        _cap_finite = (station.capacity is not None
                                       and np.isfinite(station.capacity)
                                       and station.capacity >= 0)
                        _ccap = getattr(station, '_class_capacity', None) or {}
                        _jc = self._classes[r] if r < len(self._classes) else None
                        _ccap_v = _ccap.get(_jc, None)
                        # class_cap==0 is the class-absent sentinel, not a real buffer; mirrors MATLAB refreshCapacity classCapIsFinite fix.
                        _ccap_finite = (_ccap_v is not None and np.isfinite(_ccap_v)
                                        and _ccap_v > 0)
                        if _cap_finite or _ccap_finite:
                            raise RuntimeError(
                                "Station '%s' declares setDropRule(WAITQ) for the open class '%s' "
                                "at a finite capacity. LINE does not implement waiting-room blocking "
                                "for an open arrival at a plain finite buffer: no solver honours this "
                                "combination (SolverCTMC drops the arrival, SolverJMT blocks at the "
                                "Source and ignores the capacity). Use DropStrategy.DROP for a loss "
                                "station (M/M/1/K), or one of the blocking policies LINE implements "
                                "(DropStrategy.BAS, DropStrategy.BBS, DropStrategy.RSRD) for blocking "
                                "between stations."
                                % (station.getName(), _jc.getName() if _jc is not None else str(r + 1)))

                if station_drop_rule is None:
                    # No explicit rule - default based on capacity and class type,
                    # as in MATLAB refreshCapacity.m.
                    _njobs_r = float(np.asarray(self._sn.njobs, dtype=float).ravel()[r]) \
                        if getattr(self._sn, 'njobs', None) is not None else np.inf
                    if not (np.isfinite(station.capacity) and station.capacity >= 0):
                        # No buffer: the rule is never consulted.
                        droprule[ist, r] = SnDropStrategy.WAITQ
                    elif not np.isinf(_njobs_r):
                        # finite buffer + closed class keeps WAITQ (DROP would destroy jobs; JMT and CTMC disagree on losing them).
                        droprule[ist, r] = SnDropStrategy.WAITQ
                    else:
                        # Real finite buffer and an OPEN class: LINE's reference
                        # (CTMC) loses the arrival, so DROP.
                        droprule[ist, r] = SnDropStrategy.DROP
                else:
                    # DropStrategy ids agree with lang/base.py; an unmappable rule must raise, not silently become DROP.
                    droprule[ist, r] = int(station_drop_rule)

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

        # node._prob_routing merged into _routing_matrix; skipped when _original_routes already reflects post-CS-insertion routes.
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

        # closed-class jobs start ONLY at their reference station (not every Delay node with nonzero population).
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

        # default RAND routing reaches all reachable nodes, open and closed alike; mirrors MATLAB getRoutingMatrix.m:119-136.
        if hasattr(self._sn, 'connmatrix') and self._sn.connmatrix is not None:
            connmatrix = self._sn.connmatrix
            for ind in range(nnodes):
                # If this node has any incoming connection, add it to has_incoming_route
                # for all classes (both closed and open)
                if np.sum(connmatrix[:, ind]) > 0:
                    for k in range(nclasses):
                        has_incoming_route.add((ind, k))

        # default RAND routing applied for classes at connected nodes lacking explicit routing; mirrors MATLAB getRoutingMatrix.m:119-136.
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

        # hit/miss classes only exist at and downstream of their Cache node.
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

                # skip classes with no incoming route to this node, except Source (all open classes) and ClassSwitch (all classes).
                if (ind, k) not in has_incoming_route and not is_source and not is_classswitch:
                    continue

                # skip hit/miss classes at non-Cache nodes, except via ClassSwitch/Router pass-through.
                if k in hit_miss_class_indices and not is_cache:
                    # Check if this is a ClassSwitch or Router node that handles this class
                    from .nodes import ClassSwitch, Router
                    if not isinstance(node, (ClassSwitch, Router)):
                        continue

                # A DISABLED class still leaves the node uniformly over its
                # connected edges, sink included: MATLAB getRoutingMatrix.m
                # gives the DISABLED case its own arm writing 1/sum(conn(ind,:))
                # ("we set this to be non-zero as otherwise the classes that do
                # not visit a classswitch are misconfigured in JMT"), which is
                # what fills the identity pass-through rows of an auto-added
                # CS_i_to_j node. Only a SelfLoopingClass and the Signal family
                # get a zero row, and that exception is what keeps a REPLY
                # signal off the nodes it never visits.
                jobclass = self._classes[k]
                strat = getattr(node, '_routing_strategies', {}).get(jobclass)
                if strat is not None:
                    strat_value = strat.value if hasattr(strat, 'value') else int(strat)
                    if strat_value == RoutingStrategy.DISABLED.value:
                        from .classes import (SelfLoopingClass as _SelfLooping,
                                              Signal as _Signal,
                                              OpenSignal as _OpenSignal,
                                              ClosedSignal as _ClosedSignal)
                        if not isinstance(jobclass, (_SelfLooping, _Signal,
                                                     _OpenSignal, _ClosedSignal)):
                            num_connections = np.sum(self._sn.connmatrix[ind, :])
                            if num_connections > 0:
                                for jnd in range(nnodes):
                                    if self._sn.connmatrix[ind, jnd] > 0:
                                        self._sn.rtnodes[ind * nclasses + k,
                                                         jnd * nclasses + k] = 1.0 / num_connections
                        continue

                # For nodes/classes without explicit routing, add default RAND routing
                # unless node has WRROBIN with explicit weights

                # Check if this node/class uses WRROBIN with explicit weights
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

        # Fork branch probabilities stay 1.0 each (not 1/fanout); row-sum >1 is normalized later in sn_refresh_visits.

        # Apply ClassSwitch matrices to rtnodes
        # This matches MATLAB's getRoutingMatrix behavior for StatelessClassSwitcher
        from .nodes import ClassSwitch, Cache

        # Cache class transformation uses a placeholder 0.5/0.5 hit/miss split, refined during analysis.
        for ind, node in enumerate(self._nodes):
            if isinstance(node, Cache):
                # Get hit/miss class mappings
                hit_classes = getattr(node, '_hit_class', {})  # Dict[input_class, output_class]
                miss_classes = getattr(node, '_miss_class', {})  # Dict[input_class, output_class]
                rclasses = getattr(node, '_retrieval_classes', {})  # {(item,input):retrievalClass}

                if not (hit_classes or miss_classes):
                    continue

                # Extract routing rows for this node (snapshot before transform)
                Pi = self._sn.rtnodes[ind * nclasses:(ind + 1) * nclasses, :].copy()

                def _cidx(c):
                    return c._index if hasattr(c, '_index') else self._classes.index(c)

                if rclasses:
                    # retrieval-system cache: requesting class switches to hit or a retrieval class, each routed via its own row; see _kb/09-ldes-and-cache.md.
                    for input_class, hit_class in hit_classes.items():
                        input_idx = _cidx(input_class)
                        outputs = []
                        if hit_class is not None:
                            outputs.append(_cidx(hit_class))
                        for (it, ic), oc in rclasses.items():
                            if ic is input_class and oc is not None:
                                outputs.append(_cidx(oc))
                        outputs = [o for o in outputs if o >= 0]
                        if not outputs:
                            continue
                        share = 1.0 / len(outputs)
                        self._sn.rtnodes[ind * nclasses + input_idx, :] = 0
                        for out_idx in outputs:
                            for jnd in range(nnodes):
                                out_route = Pi[out_idx, jnd * nclasses + out_idx]
                                if out_route > 1e-10:
                                    self._sn.rtnodes[ind * nclasses + input_idx,
                                                     jnd * nclasses + out_idx] += share * out_route
                    # a returning retrieval class re-enters the cache as itself; the miss is logged there, no separate routing needed.
                else:
                    # Plain cache: split the requesting class 0.5/0.5 between its
                    # hit and miss classes, following the requesting class's row.
                    for input_class, hit_class in hit_classes.items():
                        input_idx = _cidx(input_class)
                        hit_idx = _cidx(hit_class) if hit_class is not None else -1
                        miss_class = miss_classes.get(input_class)
                        miss_idx = _cidx(miss_class) if miss_class is not None else -1

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

        # open classes routed Sink->Source to close the loop for stationary computation; mirrors MATLAB getRoutingMatrix.m:325-331.
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

        # rt = stochastic complement of rtnodes, Sink->Source folded FIRST, else open rows substochastic (unbounded fluid-closing mass); mirrors MATLAB/JAR.
        from ..api.mc.dtmc import dtmc_stochcomp

        # stateful node index list mirrors MATLAB find(sn.isstateful).
        stateful_nodes_classes = []
        if self._sn.isstateful is not None:
            stateful_node_indices = np.where(self._sn.isstateful)[0]
            for ind in stateful_node_indices:
                for k in range(nclasses):
                    stateful_nodes_classes.append(int(ind) * nclasses + k)

        if len(stateful_nodes_classes) > 0:
            stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

        # Add Sink->Source routing to rtnodes, closing the open classes into a
        # loop that makes the DTMC ergodic for both rt and the visit ratios.
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
                else:
                    # All rates are zero (e.g., all arrivals disabled): use equal probabilities
                    # to avoid NaN/missing entries in Sink->Source routing
                    n_chain = len(others_in_chain)
                    for k in others_in_chain:
                        for m in others_in_chain:
                            self._sn.rtnodes[idx_sink * nclasses + k, idx_source * nclasses + m] = 1.0 / n_chain

        # rt is indexed by stateful nodes (nstateful*nclasses)^2, absorbing inserted ClassSwitch nodes.
        if len(stateful_nodes_classes) > 0:
            self._sn.rt = dtmc_stochcomp(self._sn.rtnodes, stateful_nodes_classes)

        # rt already carries the Sink->Source closure; the visit-ratio matrix field is kept only for explicit callers.
        self._sn.rt_visits = self._sn.rt

        # Connection matrix was already built earlier in this method

        # routing strategy matrix built late, defaulting to RAND before sn.routing is available.
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
            # explicit set_prob_routing sets strategy PROB unless an explicit non-PROB strategy is set; mirrors MATLAB setProbRouting.
            if hasattr(node, '_prob_routing') and node._prob_routing:
                for jobclass in node._prob_routing.keys():
                    class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                    # Don't override an explicit strategy; RAND was demoted here, and JMT writes Random for it and Empirical for PROB, which draw differently.
                    if hasattr(node, '_routing_strategies') and jobclass in node._routing_strategies:
                        explicit = node._routing_strategies[jobclass]
                        ev = explicit.value if hasattr(explicit, 'value') else int(explicit)
                        if ev != int(RoutingStrategy.PROB):
                            continue  # keep the explicit strategy
                    self._sn.routing[i, class_idx] = int(RoutingStrategy.PROB)

            # ClassSwitch nodes force PROB routing so the switch matrix's probabilities apply (RAND would be uniform and wrong).
            from .nodes import ClassSwitch
            if isinstance(node, ClassSwitch):
                for k in range(nclasses):
                    self._sn.routing[i, k] = int(RoutingStrategy.PROB)

        # Build routing weights dictionary for WRROBIN routing
        # Dict mapping (node_idx, class_idx) -> {dest_node_idx: weight}
        self._sn.routingweights = {}
        for i, node in enumerate(self._nodes):
            if hasattr(node, '_routing_weights') and node._routing_weights:
                for jobclass, dest_weights in node._routing_weights.items():
                    class_idx = jobclass._index if hasattr(jobclass, '_index') else self._classes.index(jobclass)
                    key = (i, class_idx)
                    self._sn.routingweights[key] = {}
                    for dest_node, weight in dest_weights.items():
                        # Get destination node index
                        dest_idx = dest_node._index if hasattr(dest_node, '_index') else self._nodes.index(dest_node)
                        self._sn.routingweights[key][dest_idx] = weight

    def _holdsReplyFor(self, sn, ind: int, r: int) -> bool:
        """Whether node `ind` holds a server across a synchronous call whose
        reply class is `r`, i.e. `r` returns there to RELEASE a server rather
        than to be served by one. Both indices are 0-based.

        `sn.syncreply` is indexed by the CALLING class and holds the 0-based
        reply class, -1 where no reply is expected; `sn.replyblock` marks the
        (node, calling class) pairs that hold a server across the call.
        """
        block = getattr(sn, 'replyblock', None)
        reply = getattr(sn, 'syncreply', None)
        if block is None or reply is None:
            return False
        block = np.asarray(block)
        reply = np.asarray(reply).ravel()
        if block.ndim != 2 or ind >= block.shape[0]:
            return False
        for k in range(min(reply.size, block.shape[1])):
            if int(reply[k]) == r and block[ind, k] != 0:
                return True
        return False

    def _check_service_reachable(self) -> None:
        """Refuse a class that is ROUTED TO a station which cannot serve it.

        `sanitize` disables the OUTGOING routing of a class a station cannot
        serve, which is what keeps it out of that station's visit ratios -- but
        nothing stopped the class being routed IN, and a class that arrives
        where it cannot be served is a flow sink: it enters and never leaves.
        The existing guard could not see this because it asks whether the
        station serves ANY class, not whether it serves the classes that reach
        it.

        What the three solvers did with such a model, on a closed cycle
        D <-> Q whose class C2 has no service at Q, is why this raises rather
        than warns -- one model, three different wrong answers, none flagged:

            MVA   Q/C2 ArvR 1, Tput 0        (flow not conserved)
            CTMC  drops class C2 entirely
            SSA   D/C2 QLen 2e-06, Q rows absent

        An absent slot and an explicit `Disabled()` are treated ALIKE. The
        distinction is real at the model level -- `Disabled()` declares "this
        class does not come here" -- but it does not separate the broken models
        from the sound ones: `sanitize` fills an absent slot with `Disabled()`
        and both then produce the identical struct and the identical wrong
        numbers. What separates them is whether anything actually routes the
        class in, which is what this tests. A class marked `Disabled()` at a
        station it never visits -- the normal marker in a class-switching model,
        where each queue serves exactly one class -- has no incoming flow there
        and is untouched.

        READS `sn.rtnodes`, WALKED FORWARD FROM THE FEED POINTS -- not
        `sn.nodevisits`, which this guard read until 2026-09-02 and which
        94d5570f3 had made blind to the very case it exists for. That commit
        extended the `served` mask of sn_refresh_visits from the station chain
        to the NODE chain, and it had to: on a materialised LQN replica the
        unserved states close into a whole spurious cycle. But the mask zeroes
        exactly the (station, class) cell a flow sink shows up in. On
        Source -> Q -> Sink with class B unservable at Q, B's chain went from
        Q = 1 to Q = 0 and the guard fell silent, while the Sink still read 1 --
        flow arriving downstream of a node it never visited. A MASKED VISIT
        VECTOR CANNOT ANSWER THIS QUESTION, because the mask IS the answer being
        looked for. Do not route this guard back through `nodevisits` or
        `visits`; both carry that mask.

        `rtnodes` on its own over-approximates -- it says where a class WOULD go
        if one existed -- and the WALK is what removes the slack. It starts only
        at (Source, class) pairs whose arrival is not Disabled, and at the
        reference station of each closed class with a positive population, so
        the Disabled-arrival row a class-switching Source carries is never
        entered. Three rules keep it honest, each answering a case the old visit
        solve got right for free:

        - A SINK IS ABSORBING. `rtnodes` wraps the Sink back to the Source to
          close the kernel; following that wrap re-enters every Source row,
          including the Disabled ones the seeding just excluded.
        - A SOURCE IS EXPANDED ONLY AS A SEED, for that same reason.
        - AN UNSERVABLE (station, class) IS REACHED BUT NOT EXPANDED. Nothing
          leaves it -- that is the whole complaint -- so nothing downstream of
          it is evidence of anything.

        This SUBSUMES the fed-chain precondition the guard used to carry
        separately: a chain no job can enter has no seed, so its rows are never
        walked at all. That is strictly finer than the per-chain test it
        replaces, which admitted every class of a chain any one of whose classes
        was fed.

        A SYNCHRONOUS REPLY IS NOT SERVED BY THE STATION IT RETURNS TO: it
        releases the server that station held across the call, which is the
        whole content of set_sync_reply. Its Disabled service there is the
        marker of the feature, not a flow sink, so the (station, reply class)
        pairs `sn.replyblock` marks are exempt.
        """
        if not getattr(self, '_do_checks', True):
            return
        sn = self._sn
        reached = self._reachedNodeClasses(sn)
        if reached is None or reached.size == 0:
            return
        from ..distributions import Disabled
        from .nodes import Cache, Delay, Fork, Join, Place, Queue, Transition
        # Same exemption as the sanitize checks: a station of a cache, Petri-net
        # or fork-join model legitimately carries no per-class service.
        if any(isinstance(nd, (Cache, Place, Transition, Fork, Join))
               for nd in self._nodes):
            return
        R = len(self._classes)
        if R == 0:
            return
        for ind, node in enumerate(self._nodes):
            if not isinstance(node, (Queue, Delay)):
                continue
            if getattr(node, '_hetero_service', None):
                continue
            store = getattr(node, '_service_process', None)
            if store is None:
                continue
            if ind >= reached.shape[0]:
                continue
            for r, k in enumerate(self._classes):
                svc = store.get(k)
                if svc is not None and not isinstance(svc, Disabled):
                    continue
                if r >= reached.shape[1] or not reached[ind, r]:
                    continue
                # A SYNCHRONOUS REPLY is not served by the station it returns
                # to: it releases the server that station held across the call,
                # which is the whole content of set_sync_reply.
                if self._holdsReplyFor(sn, ind, r):
                    continue
                raise ValueError(
                    "[{0}] {1} '{2}' has no service configured for job class "
                    "'{3}', but the class is routed to it. Jobs would arrive "
                    "and never leave. Call setService() for that class, or "
                    "route it elsewhere.".format(
                        self.name,
                        'Delay' if isinstance(node, Delay) else 'Queue',
                        node.getName(), k.getName()))

    def _servesClass(self, sn, ind: int, node, r: int, k) -> bool:
        """Whether a job of class `r` can LEAVE node `ind` again.

        True for anything that is not a service station, and for a station that
        serves `r`, declares heterogeneous server types (its per-class slot is
        empty by construction) or holds a server across a synchronous call whose
        reply class is `r`. False only for the flow sink itself, which is what
        stops the walk in `_reachedNodeClasses` -- the same three exemptions the
        guard applies, kept in one place so the walk and the verdict cannot
        drift apart.
        """
        from ..distributions import Disabled
        from .nodes import Delay, Queue
        if not isinstance(node, (Queue, Delay)):
            return True
        if getattr(node, '_hetero_service', None):
            return True
        store = getattr(node, '_service_process', None)
        if store is None:
            return True
        svc = store.get(k)
        if svc is not None and not isinstance(svc, Disabled):
            return True
        return self._holdsReplyFor(sn, ind, r)

    def _reachedNodeClasses(self, sn):
        """The (node, class) pairs a job can actually ARRIVE at, as an (N, R) bool array.

        A forward walk of `sn.rtnodes` from the feed points. See
        `_check_service_reachable` for why the evidence is the UNMASKED routing
        kernel rather than `sn.nodevisits`, and for the three rules -- absorbing
        Sink, Source expanded only as a seed, unservable pair reached but not
        expanded -- that keep the walk from over-approximating.

        Returns None when `rtnodes` is unreadable or not the expected
        (N*R, N*R), which leaves the caller checking nothing, as before.
        """
        from ..distributions import Disabled
        from .classes import ClosedClass
        from .nodes import Sink, Source
        rt = getattr(sn, 'rtnodes', None)
        if rt is None:
            return None
        try:
            rt = np.asarray(rt, dtype=float)
        except (TypeError, ValueError):
            return None
        N = len(self._nodes)
        R = len(self._classes)
        if N == 0 or R == 0 or rt.ndim != 2 or rt.shape[0] < N * R or rt.shape[1] < N * R:
            return None
        reached = np.zeros((N, R), dtype=bool)
        seeds = set()
        for ind, node in enumerate(self._nodes):
            if not isinstance(node, Source):
                continue
            store = getattr(node, '_arrival_process', {}) or {}
            for r, k in enumerate(self._classes):
                dist = store.get(k)
                if dist is not None and not isinstance(dist, Disabled):
                    seeds.add((ind, r))
        for r, k in enumerate(self._classes):
            if not isinstance(k, ClosedClass):
                continue
            try:
                if float(k.getNumberOfJobs()) <= 0:
                    continue
            except (TypeError, ValueError):
                continue
            ref = getattr(k, '_refstat', None)
            if ref is None:
                continue
            for ind, node in enumerate(self._nodes):
                if node is ref:
                    seeds.add((ind, r))
                    break
        stack = []
        for ind, r in seeds:
            if not reached[ind, r]:
                reached[ind, r] = True
                stack.append((ind, r))
        while stack:
            ind, r = stack.pop()
            node = self._nodes[ind]
            if isinstance(node, Sink):
                continue
            if isinstance(node, Source) and (ind, r) not in seeds:
                continue
            if not self._servesClass(sn, ind, node, r, self._classes[r]):
                continue
            row = rt[ind * R + r]
            for t in np.nonzero(row[:N * R] > GlobalConstants.Zero)[0]:
                j, sIdx = divmod(int(t), R)
                if not reached[j, sIdx]:
                    reached[j, sIdx] = True
                    stack.append((j, sIdx))
        return reached

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

        # csmask (class transition graph) built from rtnodes, before absorption, so non-station switches (ClassSwitch, Cache) are captured.
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

        # Persist the class-switching mask: csmask[r, s] is True if class r can
        # switch into class s at some node (mirrors MATLAB sn.csmask).
        self._sn.csmask = class_connected

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

        # chains formed by union over csmask (class-switching reachability), NOT by shared reference station; mirrors MATLAB getRoutingMatrix.m:279-291.
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

        # All classes in a chain must share one reference station
        # (matches the hard error in MATLAB refreshRoutingMatrix.m)
        refstat = self._sn.refstat.flatten() if self._sn.refstat is not None else np.zeros(nclasses, dtype=int)
        for chain_id in range(nchains):
            classes_in_chain = inchain[chain_id]
            if len(classes_in_chain) > 0:
                first_refstat = int(refstat[classes_in_chain[0]]) if classes_in_chain[0] < len(refstat) else 0
                for k in classes_in_chain:
                    if k < len(refstat) and int(refstat[k]) != first_refstat:
                        raise ValueError(
                            f"Classes within chain {chain_id} "
                            f"(classes: {[int(x) for x in classes_in_chain]}) "
                            f"have different reference stations.")

        # Compute visit ratios for each chain
        self._sn.visits = {}

        # visits computed from rt_visits (stochastic complement w/ Sink->Source), matching MATLAB sn_refresh_visits.m Pchain.
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
                                            # cache hit/miss class switching, folded into rt, is redistributed proportionally over rt[src,:] for classes h and m.
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

            # sparsify once; the per-chain routing matrices are slices of this
            from ..api.sn.transforms import _csr_normalize_rows, _csr_threshold_mask
            rt_sparse = sp.csr_matrix(rt)

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

                # chain submatrix indexed by stateful nodes (M=nstateful), not by all nodes; mirrors MATLAB Pchain = rt(cols,cols).
                nstateful = self._sn.nstateful

                # Build indices in rt space (stateful_idx * nclasses + class_idx)
                indices = []
                for isf in range(nstateful):
                    for k in classes_in_chain:
                        indices.append(isf * nclasses + k)

                if len(indices) > 0:
                    # the per-chain routing matrix is kept in sparse storage
                    P_chain = sp.csr_matrix(rt_sparse[indices, :][:, indices])

                    # Use DTMC solver for both open and closed chains
                    # This matches MATLAB's sn_refresh_visits.m which uses dtmc_solve for all chains
                    from ..api.mc.dtmc import dtmc_solve, dtmc_solve_reducible
                    from .base import NodeType

                    n = len(indices)

                    # the routing matrix carries a JMT-oriented uniform fill on DISABLED (node,class)
                    # pairs; a class with no service at a station cannot be there -- see _kb/06-solver-catalog.md
                    if self._sn.rates is not None:
                        from .base import NodeType as _NTv
                        from ..api.sn.transforms import _sn_has_server_types
                        _served = np.ones(P_chain.shape[0])
                        for isf in range(nstateful):
                            _sti = int(self._sn.statefulToStation[isf]) if self._sn.statefulToStation is not None else -1
                            if _sti < 0 or _sti >= self._sn.nstations:
                                continue
                            # a Place, and a station declaring server types, carry NaN station rates by construction
                            if self._sn.nodetype[int(self._sn.stationToNode[_sti])] in (_NTv.PLACE, _NTv.TRANSITION):
                                continue
                            if _sn_has_server_types(self._sn, _sti):
                                continue
                            for _ik, _k in enumerate(classes_in_chain):
                                if not np.isnan(self._sn.rates[_sti, int(_k)]):
                                    continue
                                _idx = isf * n_chain_classes + _ik
                                if _idx < _served.shape[0]:
                                    _served[_idx] = 0.0
                        if _served.size and not _served.all():
                            _D = sp.diags(_served)
                            P_chain = sp.csr_matrix(_D @ P_chain @ _D)
                            P_chain.eliminate_zeros()

                    row_sums = np.asarray(P_chain.sum(axis=1)).ravel()
                    visited = row_sums > 0

                    # Check for fork nodes
                    has_fork = False
                    if hasattr(self._sn, 'nodetype') and self._sn.nodetype is not None:
                        fork_type = NodeType.FORK if hasattr(NodeType, 'FORK') else getattr(NodeType, 'Fork', None)
                        if fork_type is not None:
                            # Use Python any() for enum comparison (np.any doesn't work correctly with enum lists)
                            has_fork = any(nt == fork_type for nt in self._sn.nodetype)

                    # open chains with zero total arrival rate get zero visits (mirrors MATLAB's 0/0->NaN->zero propagation through the DTMC solve).
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

                    # open fork networks use the same DTMC path as closed forks now that rt_visits folds Sink->Source routing.

                    if np.sum(visited) > 0:
                        _vidx = np.where(visited)[0]
                        P_visited = P_chain[_vidx, :][:, _vidx]

                        # Fork rows (sum>1) normalized before dtmc_solve_reducible, with fanout correction applied after; mirrors MATLAB sn_refresh_visits.m:88-98.
                        row_sums_visited = np.ones(P_visited.shape[0])
                        if has_fork:
                            row_sums_visited = _csr_normalize_rows(P_visited, 1e-10)

                        # dtmc_solve is primary, dtmc_solve_reducible the fallback for reducible chains; see _kb/11-conventions-and-gotchas.md DTMC solver order.
                        from scipy.sparse.csgraph import connected_components
                        from scipy.sparse import csc_matrix
                        n_components, _ = connected_components(
                            csc_matrix(_csr_threshold_mask(P_visited, 1e-10)), directed=True, connection='strong', return_labels=True
                        )
                        is_reducible = n_components > 1

                        # dtmc_solve PRIMARY, dtmc_solve_reducible fallback -- order is load-bearing, never change it; see _kb/11-conventions-and-gotchas.md.
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

                        # SPN-based fork correction replaces the transitive-closure correction: population-preserving SPN analysis gives uniform visit ratios on fork-join.
                        if has_fork and np.any(row_sums_visited > 1 + 1e-10):
                            for idx in range(len(pi)):
                                if pi[idx] > 1e-10:
                                    pi[idx] = 1
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

                    # normalize by the SUM of visits at the chain's reference station over all its classes; mirrors MATLAB sn_refresh_visits.m:118-121.
                    ref_class = classes_in_chain[0]
                    ref_stat = int(self._sn.refstat[ref_class])  # station index
                    ref_sf = int(self._sn.stationToStateful[ref_stat])  # stateful index
                    norm_sum = np.sum(visits_chain[ref_sf, classes_in_chain])
                    if norm_sum > 1e-10:
                        visits_chain = visits_chain / norm_sum
                    elif is_reducible and has_fork and is_open_chain:
                        # open fork networks with absorbing states: each forked branch gets visits=1/fanout directly from routing structure.
                        visits_chain = self._compute_fork_visits(
                            nstateful, nclasses, classes_in_chain, ref_sf
                        )
                    # reference station with zero visits stays exact zero, not uniform; see _kb/04-networkstruct.md Python native network.py section.
                    visits_chain = np.abs(visits_chain)
                    visits_chain[visits_chain < GlobalConstants.Zero] = 0.0

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
            from ..api.sn.transforms import _csr_fill_nan_rows, _csr_normalize_rows

            # sparsify once; the per-chain node routing matrices are slices of this
            rtnodes_sparse = sp.csr_matrix(self._sn.rtnodes)

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
                    # the per-chain node routing matrix is kept in sparse storage
                    nodes_Pchain = sp.csr_matrix(rtnodes_sparse[nodes_cols, :][:, nodes_cols])

                    # Handle NaN values in routing matrix
                    _csr_fill_nan_rows(nodes_Pchain)

                    # THE SAME DISABLED-PAIR MASK THE STATION BLOCK APPLIES ABOVE,
                    # and it matters more here: at station level a (station,class)
                    # the class cannot be served at is a dead end, while the node
                    # kernel keeps the class-switch nodes between the stations, so
                    # the disabled states close into a whole spurious CYCLE. A
                    # materialised LQN replica is exactly that -- replica 2's
                    # stations still carry replica 1's classes in rtnodes -- and
                    # dtmc_solve_reducible then splits the mass between the real
                    # chain and the phantom one, giving every node of replica 2 a
                    # visit in replica 1's classes.
                    if self._sn.rates is not None:
                        from .base import NodeType as _NTn
                        from ..api.sn.transforms import _sn_has_server_types as _has_st
                        _nserved = np.ones(nodes_Pchain.shape[0])
                        for _ind in range(nnodes):
                            _sti = int(self._sn.nodeToStation[_ind]) if self._sn.nodeToStation is not None else -1
                            if _sti < 0 or _sti >= self._sn.nstations:
                                continue
                            # a Place, and a station declaring server types, carry NaN station rates by construction
                            if self._sn.nodetype[_ind] in (_NTn.PLACE, _NTn.TRANSITION):
                                continue
                            if _has_st(self._sn, _sti):
                                continue
                            for _ik, _k in enumerate(classes_in_chain):
                                if not np.isnan(self._sn.rates[_sti, int(_k)]):
                                    continue
                                _idx = _ind * n_chain_classes + _ik
                                if _idx < _nserved.shape[0]:
                                    _nserved[_idx] = 0.0
                        if _nserved.size and not _nserved.all():
                            _Dn = sp.diags(_nserved)
                            nodes_Pchain = sp.csr_matrix(_Dn @ nodes_Pchain @ _Dn)
                            nodes_Pchain.eliminate_zeros()

                    # Find visited nodes
                    nodes_visited = np.asarray(nodes_Pchain.sum(axis=1)).ravel() > 0

                    if np.sum(nodes_visited) > 0:
                        _nvidx = np.where(nodes_visited)[0]
                        nodes_Pchain_visited = nodes_Pchain[_nvidx, :][:, _nvidx]

                        # Normalize rows, recording original row sums for fork correction
                        # (matches MATLAB sn_refresh_visits.m lines 198-209)
                        nodes_row_sums_visited = _csr_normalize_rows(nodes_Pchain_visited, 1e-10)

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

                        # SPN-based fork correction for node visits: stations/Fork get visit=1,
                        # Join nodes get visit = number of direct predecessors.
                        if has_fork and np.any(nodes_row_sums_visited > 1 + 1e-10):
                            for idx in range(len(nodes_alpha)):
                                if nodes_alpha[idx] > 1e-10:
                                    nd = idx // n_chain_classes
                                    join_type = NodeType.JOIN if hasattr(NodeType, 'JOIN') else getattr(NodeType, 'Join', None)
                                    if nd < len(self._sn.nodetype) and join_type is not None and self._sn.nodetype[nd] == join_type:
                                        r = classes_in_chain[idx % n_chain_classes]
                                        col = nd * nclasses + r
                                        if col < self._sn.rtnodes.shape[1]:
                                            n_sources = int(np.sum(self._sn.rtnodes[:, col] > 1e-10))
                                            nodes_alpha[idx] = n_sources
                                        else:
                                            nodes_alpha[idx] = 1
                                    else:
                                        nodes_alpha[idx] = 1
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

                    # normalize by reference-station visits via statefulToNode(refstat(...)); mirrors MATLAB sn_refresh_visits.m:220 exactly for FJ aux-chain parity.
                    ref_class = classes_in_chain[0]
                    ref_stat = int(self._sn.refstat[ref_class])
                    if ref_stat < len(self._sn.statefulToNode):
                        ref_node = self._sn.statefulToNode[ref_stat]
                    else:
                        ref_node = 0
                    node_norm_sum = np.sum(nodevisits_chain[ref_node, classes_in_chain])
                    if node_norm_sum > 1e-10:
                        nodevisits_chain = nodevisits_chain / node_norm_sum

                    nodevisits_chain[nodevisits_chain < 0] = 0
                    nodevisits_chain[np.isnan(nodevisits_chain)] = 0

                    self._sn.nodevisits[chain_id] = nodevisits_chain

        # refclass=-1 means no reference class (python 0-indexed vs MATLAB's refclass=0); used as the WN divisor in sn_get_residt_from_respt.
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
                    # per-mode SPN firing process: (D0,D1) Markovian phase-type, phases NaN for non-Markovian pending sn_nonmarkov_toph, pie/procid also recorded.
                    firingproc: List[Any] = field(default_factory=list)
                    firingphases: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
                    firingpie: List[Any] = field(default_factory=list)
                    firingprocid: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
                    # Marking-dependent firing-rate multiplier per mode; None == unit.
                    firingdep: List[Any] = field(default_factory=list)

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

                # Marking-dependent firing-rate multiplier per mode (None == unit)
                frm = getattr(node, '_firing_rate_dependence', None)
                param.firingdep = [
                    (frm[mode_idx] if frm is not None and mode_idx < len(frm) else None)
                    for mode_idx in range(nmodes)
                ]

                # phase-type firing process per mode; mirrors MATLAB refreshPetriNetNodes.m:47-63 and JAR TransitionNodeParam.
                from ..distributions.base import Markovian as _Markovian
                from ..distributions.base import ContinuousDistribution as _ContDist
                from ..constants import ProcessType as _ProcessType
                param.firingproc = [None] * nmodes
                param.firingpie = [None] * nmodes
                param.firingphases = np.full(nmodes, np.nan, dtype=float)
                param.firingprocid = np.full(nmodes, -1, dtype=int)
                for mode_idx in range(nmodes):
                    dist = param.distributions[mode_idx] if mode_idx < len(param.distributions) else None
                    if dist is None:
                        continue
                    proc_id = _ProcessType.fromString(dist.__class__.__name__)
                    if proc_id is not None:
                        param.firingprocid[mode_idx] = proc_id.value
                    if isinstance(dist, _Markovian):
                        D0 = np.atleast_2d(np.asarray(dist.getD0(), dtype=float))
                        D1 = np.atleast_2d(np.asarray(dist.getD1(), dtype=float))
                        param.firingproc[mode_idx] = (D0, D1)
                        try:
                            pie = np.atleast_1d(np.asarray(dist.getInitProb(), dtype=float))
                        except (NotImplementedError, AttributeError):
                            pie = np.zeros(D0.shape[0], dtype=float)
                            if pie.size > 0:
                                pie[0] = 1.0
                        param.firingpie[mode_idx] = pie
                        param.firingphases[mode_idx] = float(D0.shape[0])
                    elif isinstance(dist, _ContDist):
                        # Non-Markovian: leave firingproc empty and firingphases NaN;
                        # sn_nonmarkov_toph will fill these in for CTMC analysis.
                        param.firingproc[mode_idx] = None
                        param.firingpie[mode_idx] = None
                        param.firingphases[mode_idx] = np.nan

                nodeparam[node_idx] = param

        # Handle Queue nodes with polling and switchover
        from .nodes import Queue
        from ..api.sn import SchedStrategy

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Queue):
                # setup/delay-off times kept as distribution objects, matching switchoverTime storage; mirrors MATLAB refreshLocalVars.m.
                if node.is_delay_off_enabled():
                    for r, jobclass in enumerate(self._classes):
                        setup = node.get_setup_time(jobclass)
                        delayoff = node.get_delay_off_time(jobclass)
                        if setup is None or delayoff is None:
                            continue
                        if node_idx not in nodeparam:
                            nodeparam[node_idx] = {}
                        if r not in nodeparam[node_idx]:
                            nodeparam[node_idx][r] = {}
                        nodeparam[node_idx][r]['setupTime'] = setup
                        nodeparam[node_idx][r]['delayoffTime'] = delayoff

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

        # Handle pass-and-swap (PAS) queues: store the class compatibility/swap
        # graph and the total service rate function mu(c) at the node level.
        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Queue) and getattr(node.get_sched_strategy(), 'name', None) in ('PAS', 'OI'):
                if node_idx not in nodeparam or not isinstance(nodeparam[node_idx], dict):
                    nodeparam[node_idx] = {}
                if getattr(node.get_sched_strategy(), 'name', None) == 'OI':
                    # Order-independent: swap graph is always zero (empty), so
                    # class order is preserved on completion (plain OI).
                    swap_graph = np.zeros((nclasses, nclasses))
                else:
                    swap_graph = node.get_swap_graph()
                    if swap_graph is None:
                        # PAS default: complete compatibility graph (no self-loops).
                        swap_graph = np.ones((nclasses, nclasses)) - np.eye(nclasses)
                nodeparam[node_idx]['swapGraph'] = np.asarray(swap_graph, dtype=float)
                nodeparam[node_idx]['svcRateFun'] = node.get_service_rate_function()

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
                    actualdelayedhitprob: Optional[np.ndarray] = None
                    actualhitproblist: Optional[np.ndarray] = None
                    actualitemprob: Optional[np.ndarray] = None
                    actualresidt: Optional[np.ndarray] = None
                    actuallistcost: Optional[np.ndarray] = None
                    # per-item storage costs and per-list cost caps (ton21cache Sec. IX)
                    itemsize: Optional[np.ndarray] = None
                    costcap: Optional[np.ndarray] = None
                    costcapglobal: bool = False
                    replacestrat: Any = None  # Renamed from 'replacement' to match MATLAB
                    qlru: float = 1.0  # q-LRU admission probability on a miss
                    accost: Optional[np.ndarray] = None
                    pread: List[Any] = field(default_factory=list)  # PMF values per class
                    # retrieval system
                    total_cache_capacity: int = 0
                    retrieval_system_capacity: int = 0
                    # [nitems x nclasses] array of 0-based class indices, -1 where undefined
                    retrieval_classes: np.ndarray = field(default_factory=lambda: np.array([]))
                    # set of 0-based retrieval class indices
                    retrieval_class_indices: set = field(default_factory=set)
                    # 0-based arrival class -> list of 0-based retrieval-system queue node indices
                    retrieval_system_queue_indices: dict = field(default_factory=dict)

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
                param.itemsize = getattr(node, '_item_size', None)
                param.costcap = getattr(node, '_cost_cap', None)
                param.costcapglobal = bool(getattr(node, '_cost_cap_global', False))
                param.replacestrat = node._replacement_strategy if hasattr(node, '_replacement_strategy') else None
                param.qlru = float(getattr(node, '_admission_prob', 1.0))

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

                # retrieval system: cache node's (item,class)->class map converted to a [nitems x nclasses] 0-based class-index array (-1 undefined).
                param.total_cache_capacity = getattr(node, '_total_cache_capacity', param.cap)
                param.retrieval_system_capacity = getattr(node, '_retrieval_system_capacity', 0)
                rc = np.full((max(param.nitems, 1), nclasses), -1, dtype=int)
                for (item, in_class), out_class in getattr(node, '_retrieval_classes', {}).items():
                    if in_class in self._classes and out_class in self._classes and item < rc.shape[0]:
                        rc[item, self._classes.index(in_class)] = self._classes.index(out_class)
                param.retrieval_classes = rc
                param.retrieval_class_indices = set(
                    getattr(node, '_retrieval_class_indices', set()))
                param.retrieval_system_queue_indices = dict(
                    getattr(node, '_retrieval_system_queue_indices', {}))

                # Access cost (if available). A per-item graph (set_access_graph)
                # is shared by all classes, mirroring MATLAB sanitize.m.
                if hasattr(node, '_accost') and node._accost is not None:
                    param.accost = np.array(node._accost)
                elif getattr(node, '_graph', None) is not None:
                    param.accost = np.array(
                        [[np.asarray(g, dtype=float) for g in node._graph]
                         for _ in range(nclasses)])

                # CLIMB solved as an equivalent FIFO cache with unit-capacity lists (exact identity); see _kb/09-ldes-and-cache.md.
                from .base import ReplacementStrategy as _RS
                if param.replacestrat == _RS.CLIMB:
                    Cclimb = int(np.sum(param.itemcap))
                    param.itemcap = np.ones(Cclimb, dtype=int)
                    param.replacestrat = _RS.FIFO
                    param.accost = None

                # SOLVER RESULTS ARE RE-DERIVED FROM THE NODE, not dropped.
                # A CacheParam is rebuilt from scratch on every refresh, so a
                # hit/miss split an analyzer had already computed vanished from
                # sn.nodeparam the moment anything refreshed the struct, and the
                # quantities derived from it (the visits carrying the hit and
                # miss classes) silently fell back to the 0.5/0.5 routing guess.
                # MATLAB copies them back off the node here, so the struct a
                # caller holds never disagrees with the node; mirrors
                # refreshLocalVars.m ("Store actual hit/miss/latency from
                # solver results"). The node is the durable home of the result.
                _ahp = getattr(node, '_actual_hit_prob', None)
                if _ahp is not None and np.size(_ahp) > 0:
                    param.actualhitprob = np.asarray(_ahp)
                    _amp = getattr(node, '_actual_miss_prob', None)
                    if _amp is not None and np.size(_amp) > 0:
                        param.actualmissprob = np.asarray(_amp)
                _adhp = getattr(node, '_actual_delayed_hit_prob', None)
                if _adhp is not None and np.size(_adhp) > 0:
                    param.actualdelayedhitprob = np.asarray(_adhp)
                _art = getattr(node, '_actual_residt', None)
                if _art is not None and np.size(_art) > 0:
                    param.actualresidt = np.asarray(_art)

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

        # Handle Join nodes: their strategy and quorum lived only on the node
        # object, so every consumer that read sn (sn_fj_validate, the exact
        # tag construction, the JMT export) saw a Join with no strategy at all
        # and treated a PARTIAL one as standard. MATLAB's refreshLocalVars puts
        # both in nodeparam; this is that twin.
        from .nodes import Join
        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Join):
                K = len(self._classes)
                join_strategy = {}
                fan_in = {}
                for r, cls in enumerate(self._classes):
                    st = node._join_strategy.get(cls, None)
                    if st is not None:
                        join_strategy[r] = st
                    req = node._required.get(cls, None)
                    if req is not None:
                        fan_in[r] = req
                entry = nodeparam.get(node_idx)
                if not isinstance(entry, dict):
                    entry = {}
                    nodeparam[node_idx] = entry
                entry['joinStrategy'] = join_strategy
                entry['fanIn'] = fan_in

        # Handle Fork nodes with tasks per link (fanOut)
        from .nodes import Fork
        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Fork):
                tasks_per_link = node.get_tasks_per_link()
                if tasks_per_link is not None:
                    tasks_per_link = np.atleast_1d(tasks_per_link)
                    # float, not int: MATLAB stores the double and lets
                    # sn_fj_validate refuse a non-integer fanout by name. An
                    # int() here turned 1.5 into 1 and answered instead.
                    fanout_val = float(tasks_per_link.flat[0]) if tasks_per_link.size > 0 else 1.0
                else:
                    fanout_val = 1.0
                # merge, do not replace: another producer may already hold this slot
                if isinstance(nodeparam.get(node_idx), dict):
                    entry = nodeparam[node_idx]
                else:
                    entry = {}
                    nodeparam[node_idx] = entry
                entry['fanOut'] = fanout_val

                # Variable forking levels. fanOut stays the scalar every
                # existing consumer reads; the three matrices beside it are
                # (nnodes x nclasses), indexed by DESTINATION NODE rather than
                # link ordinal so a relink cannot silently permute them.
                I = len(self._nodes)
                K = len(self._classes)
                connrc = np.zeros((I, K), dtype=bool)
                if self._sn.connmatrix is not None:
                    for kdest in range(I):
                        if self._sn.connmatrix[node_idx, kdest] > 0:
                            connrc[kdest, :] = True
                fan_out_link = np.zeros((I, K))
                fan_out_prob = np.zeros((I, K))
                fan_out_dist = [[None] * K for _ in range(I)]
                fan_out_link[connrc] = fanout_val
                fan_out_prob[connrc] = 1.0

                name_to_idx = {n.get_name(): i for i, n in enumerate(self._nodes)}

                def _dests(dest_name):
                    if not dest_name:
                        return [k for k in range(I) if connrc[k, :].any()]
                    if dest_name not in name_to_idx:
                        raise ValueError(
                            'Fork override names destination "%s", which is not '
                            'a node of this model.' % dest_name)
                    return [name_to_idx[dest_name]]

                for dest, cls_idx, value in node._tasks_per_link_by_dest:
                    for kdest in _dests(dest):
                        fan_out_link[kdest, cls_idx] = value
                for dest, cls_idx, value in node._branch_prob:
                    for kdest in _dests(dest):
                        fan_out_prob[kdest, cls_idx] = value
                for dest, cls_idx, dist in node._tasks_per_link_dist:
                    for kdest in _dests(dest):
                        fan_out_dist[kdest][cls_idx] = dist
                        # the scalar slot carries the mean, so a consumer that
                        # only reads fanOutLink still sees E[tasks per link]
                        fan_out_link[kdest, cls_idx] = dist.getMean()

                entry['fanOutLink'] = fan_out_link
                entry['fanOutProb'] = fan_out_prob
                entry['fanOutDist'] = fan_out_dist
                if (node._tasks_per_link_dist or node._tasks_per_link_by_dest
                        or node._branch_prob):
                    # Keep the scalar consistent with the per-link mean, so a
                    # solver that has not been taught the matrices degrades to
                    # E[.] rather than to a value the fork never emits. The
                    # branch probability is folded in here and NOT into
                    # fanOutLink, because JMT and LDES read the two separately:
                    # fanOutLink is the count GIVEN the branch fires,
                    # fanOutProb is whether it fires at all.
                    nz = (fan_out_link * fan_out_prob)[connrc]
                    if nz.size > 0:
                        entry['fanOut'] = float(nz.mean())

        # Handle Join nodes: the per-class join rule and its quorum. fanIn is
        # the JMT numRequired of a STANDARD join (-1 = every sibling),
        # joinRequired the quorum k of a PARTIAL one: the two read the same
        # field but are written to different JMT elements, so both are carried.
        from .nodes import Join
        from .base import JoinStrategy
        nclasses = len(self._classes)
        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, Join):
                join_strategy = []
                join_required = []
                for r in range(nclasses):
                    jobclass = self._classes[r]
                    join_strategy.append(node.get_strategy(jobclass))
                    req = node.get_required(jobclass)
                    join_required.append(-1 if req is None else int(req))
                nodeparam[node_idx] = {'joinStrategy': join_strategy,
                                       'fanIn': list(join_required),
                                       'joinRequired': join_required}

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

        # Skip for models created by MMT (prevents infinite recursion:
        # mmt() → nonfjmodel.refresh_struct() → _refresh_fork_join_nodevisits() → mmt() → ...)
        if getattr(self, '_skip_fj_nodevisits', False):
            return

        # Check if there are any fork-join relationships
        if sn.fj is None or not np.any(sn.fj):
            return

        # A quorum join needs no warning here: it is declared as the JoinPartial
        # feature by get_used_lang_features, so a solver that cannot serve it
        # refuses the model by name, and the MMT fixed point that can charges the
        # k-th branch completion through fj_ordstat_exp.

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

        # one iteration per MMT auxiliary chain, with an inner loop over its auxiliary classes; mirrors MATLAB refreshStruct.m:439-465.
        from .base import NodeType
        for new_chain in range(sn.nchains, fsn.nchains):
            if new_chain not in fsn.inchain or fsn.inchain[new_chain] is None:
                continue
            chain_classes = list(fsn.inchain[new_chain])
            if len(chain_classes) == 0:
                continue

            # fjclassmap is indexed by aux_idx = fsn_class - sn.nclasses in python (vs MATLAB's positional indexing).
            any_aux_idx = None
            for cls in chain_classes:
                aux_idx_candidate = cls - sn.nclasses
                if 0 <= aux_idx_candidate < len(fjclassmap):
                    any_aux_idx = aux_idx_candidate
                    break
            if any_aux_idx is None:
                continue

            orig_fork = int(forkmap[any_aux_idx])
            orig_class = int(fjclassmap[any_aux_idx])

            # Find original chain containing orig_class
            orig_chain = None
            for c in range(sn.nchains):
                if c in sn.inchain and orig_class in sn.inchain[c]:
                    orig_chain = c
                    break
            if orig_chain is None or orig_chain not in sn.nodevisits:
                continue
            if new_chain not in fsn.nodevisits:
                continue

            # Source/Sink/Fork rows zeroed in the new chain's nodevisits; Fork may already be rewritten as Router by mmt so is often left as-is, matching MATLAB.
            Vaux = fsn.nodevisits[new_chain].copy()
            for i, nt in enumerate(fsn.nodetype):
                if nt == NodeType.SOURCE or nt == NodeType.SINK or nt == NodeType.FORK:
                    Vaux[i, :] = 0

            # Keep only columns for auxiliary classes in this new chain
            Vaux = Vaux[:, chain_classes]

            # If node counts differ, map Vaux rows by node name
            # (matches MATLAB lines 445-459)
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

            # fanOut read from the original fork's nodeparam (fanout of the fork itself is already baked into fsn via mmt).
            # nodeparam holds a DICT for a Fork, so the attribute read this used
            # to do returned None on every model and pinned fanOut_val at 1: the
            # correction silently dropped tasksPerLink. Accept both shapes.
            fanOut_val = 1
            fp = sn.nodeparam.get(orig_fork) if sn.nodeparam is not None else None
            if fp is not None:
                np_fanOut = fp.get('fanOut', None) if isinstance(fp, dict) \
                    else getattr(fp, 'fanOut', None)
                if np_fanOut is not None and np_fanOut > 0:
                    fanOut_val = np_fanOut

            # original chain's pre-mmt visits captured BEFORE the inner loop writes back, so every jaux iteration reads the same X.
            X = sn.nodevisits[orig_chain].copy()

            # Apply MMT correction: one aux class jaux -> one original class j
            orig_classes_in_chain = list(sn.inchain[orig_chain])
            for jaux in range(Vaux.shape[1]):
                if jaux >= len(orig_classes_in_chain):
                    break
                j = orig_classes_in_chain[jaux]
                self._sn.nodevisits[orig_chain][:, j] = fanOut_val * (X[:, j] + Vaux[:, jaux])

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

        # reference station and fork-branch stations both get visits=1 (fork sends full rate to each branch).

        for isf in range(nstateful):
            for k in classes_in_chain:
                if isf == ref_sf:
                    # at the reference Source, only the originating class has arrivals; switched-to classes get visits=0.
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

    def _is_declared_state_node(self, node) -> bool:
        """True when this node's state row is an INPUT rather than a default.

        A partial state is dropped whole (see `_refresh_state`), which is right
        for a station whose row `initDefault` would have written anyway. Two
        kinds of node carry a row nothing else can reconstruct:

        - a PAS station, whose ORDERED placement is a required input;
        - an SPN Place, whose INITIAL MARKING is the model. Dropping it left
          `spn_open_sevenplaces` exporting a net whose downstream places never
          receive a token: JMT reported P1 throughput 1.0111 against 2.8376 and
          P5, P6 and P7 vanished from the table entirely.
        """
        if self._is_pas_node(node):
            return True
        return type(node).__name__ == 'Place'

    def _is_pas_node(self, node) -> bool:
        """A PAS station's ordered placement is an input, not a default state."""
        sched = node.get_sched_strategy() if hasattr(node, 'get_sched_strategy') else None
        if sched is None:
            return False
        # BY NAME: a node can carry either SchedStrategy enum, and the two
        # disagree numerically (PAS is 41 in one and 37 in the other).
        return str(getattr(sched, 'name', sched)) == 'PAS'

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

        # full local rows as stored on the nodes; sn.state carries these, not the
        # per-class marginal, matching MATLAB getState (State.toMarginal exists
        # precisely because sn.state{isf} is the unreduced row)
        full_rows = {}

        # a state on a strict subset of stateful nodes is not an initialization; see _kb/11-conventions-and-gotchas.md (partial initial state).
        fully_initialized = self.has_init_state()

        for node_idx, node in enumerate(self._nodes):
            if isinstance(node, StatefulNode):
                # Get stateful index
                stateful_idx = None
                if nodeToStateful is not None and node_idx < len(nodeToStateful):
                    stateful_idx = int(nodeToStateful[node_idx])

                if stateful_idx is not None and stateful_idx >= 0 and stateful_idx < len(state):
                    # Check if this node had setState() called explicitly
                    explicitly_set = fully_initialized and getattr(node, '_state_explicitly_set', False)
                    # Get node state. A PARTIAL state is dropped whole, as
                    # MATLAB drops it: getState runs initDefault when
                    # hasInitState is false, so the default marking answers and
                    # the rows that happen to be set are not mixed into it.
                    # Copying them left Q1=2 with the reference station emptied
                    # to 0, a third reading of the same model beside MATLAB's
                    # (Think=2, Q1=0) and the C++ reader's (Think=2, Q1=2).
                    # A PAS placement survives it, exactly as initDefault.m keeps
                    # the user row (its `hasUser` branch) while defaulting every
                    # other station: the ordering is a required input there, not
                    # a default, and dropping it raises the refusal instead.
                    keep_row = fully_initialized or self._is_declared_state_node(node)
                    node_state = None
                    if keep_row:
                        node_state = node.get_state() if hasattr(node, 'get_state') else node._state if hasattr(node, '_state') else None
                    if node_state is not None:
                        node_state = np.asarray(node_state).flatten()
                        if node_state.size:
                            full_rows[stateful_idx] = node_state.copy()
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
            # closed classes w/o explicit state: all jobs at the reference station if it has capacity, else spill in index order.
            nstations_loc = len(self._stations)
            classcap_loc = getattr(self._sn, 'classcap', None)
            cap_loc = getattr(self._sn, 'cap', None)
            placed = np.zeros((nstations_loc, nclasses))
            for r in range(nclasses):
                if not explicit_state_set[r] and np.isfinite(njobs[r]) and njobs[r] > 0:
                    ref_station = int(refstat[r])
                    if ref_station >= len(stationToStateful):
                        continue
                    if classcap_loc is None or cap_loc is None:
                        placed[ref_station, r] = njobs[r]
                        continue
                    cap_flat = np.asarray(cap_loc).flatten()
                    remaining = float(njobs[r])
                    for jst in [ref_station] + [j for j in range(nstations_loc) if j != ref_station]:
                        if remaining <= 0:
                            break
                        avail = min(classcap_loc[jst, r] - placed[jst, r],
                                    cap_flat[jst] - np.sum(placed[jst, :]))
                        take = min(remaining, max(0.0, avail))
                        placed[jst, r] += take
                        remaining -= take
                    if remaining > 0:
                        # infeasible placement is caught downstream; put the rest at ref
                        placed[ref_station, r] += remaining
            for ist in range(nstations_loc):
                for r in range(nclasses):
                    if placed[ist, r] > 0:
                        stateful_idx = int(stationToStateful[ist])
                        if 0 <= stateful_idx < len(state):
                            state[stateful_idx][r] = placed[ist, r]

        # sn.state{isf} is the FULL local row (buffer order, server phases), never
        # the per-class marginal: State.toMarginal reduces it on demand, and the
        # CTMC seed lookup matches it against a row of sn.space. Storing the
        # marginal here made that lookup fail on every ordered-buffer station, which
        # silently skipped the reachability pruning and left a reducible generator
        # to be split equally across its recurrent classes. See _kb/11.
        if stationToStateful is not None:
            from ..api.state.marginal import fromMarginal
            nstations = len(self._stations) if hasattr(self, '_stations') else 0
            stationToNode = getattr(self._sn, 'stationToNode', None)
            stationToNode = np.asarray(stationToNode).flatten() if stationToNode is not None else None
            sched = getattr(self._sn, 'sched', None)
            for ist in range(nstations):
                isf = int(stationToStateful[ist])
                if isf < 0 or isf >= len(state):
                    continue
                if isf in full_rows:
                    state[isf] = full_rows[isf]
                    continue
                # No row on the node: the model was never initialized, so expand the
                # placement marginal the way initDefault would. PAS is excluded because
                # its ordering is a required input rather than a default -- expanding it
                # would fabricate the placement that picks the recurrent component, and
                # silence the refusal that must be raised instead.
                if stationToNode is None or ist >= len(stationToNode):
                    continue
                sched_i = (sched.get(ist) if hasattr(sched, 'get') else
                           (sched[ist] if sched is not None else None))
                if sched_i is not None and sched_i == SchedStrategy.PAS:
                    continue
                # Only the seed row is kept, so ask for it alone: enumerating the
                # whole local space here is binomial in the population and
                # exhausts memory before any solver gate can price the chain.
                try:
                    space_i = fromMarginal(self._sn, int(stationToNode[ist]),
                                           state[isf], top_row_only=True)
                except Exception:
                    continue
                if space_i is None:
                    continue
                space_i = np.atleast_2d(np.asarray(space_i, dtype=float))
                if space_i.size:
                    state[isf] = space_i[0, :].copy()

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

        # initial marginal: all closed jobs at reference station, spilling to others if capacity < N; Place-referenced classes keep all-at-ref (SPN tokens).
        n0 = np.zeros((nstations, nclasses))
        sn = self.get_struct()
        classcap = sn.classcap if getattr(sn, 'classcap', None) is not None \
            else np.full((nstations, nclasses), np.inf)
        cap = sn.cap if getattr(sn, 'cap', None) is not None \
            else np.full(nstations, np.inf)
        cap = np.asarray(cap).flatten()

        def _enum_val(x):
            return int(x.value) if hasattr(x, 'value') else int(x)

        def _is_place(jst):
            node_idx = int(np.asarray(sn.stationToNode).flatten()[jst])
            return _enum_val(sn.nodetype[node_idx]) == _enum_val(NodeType.PLACE)

        def _is_ext(jst):
            sched = sn.sched.get(jst) if hasattr(sn.sched, 'get') else sn.sched[jst]
            if sched is None:
                return False
            return _enum_val(sched) == _enum_val(SchedStrategy.EXT)

        for r, jobclass in enumerate(self._classes):
            # Check if this is a closed class with a reference station
            if hasattr(jobclass, '_refstat') and jobclass._refstat is not None:
                refstat = jobclass._refstat
                # Find the station index
                try:
                    ist = self._stations.index(refstat)
                except ValueError:
                    continue  # Station not found
                njobs = getattr(jobclass, '_njobs', 0)
                if not np.isfinite(njobs):
                    continue
                if _is_place(ist):
                    n0[ist, r] = njobs
                    continue
                remaining = njobs
                for jst in [ist] + [j for j in range(nstations) if j != ist]:
                    if remaining <= 0:
                        break
                    if _is_ext(jst) or _is_place(jst):
                        continue
                    avail = min(classcap[jst, r] - n0[jst, r], cap[jst] - np.sum(n0[jst, :]))
                    take = min(remaining, max(0.0, avail))
                    n0[jst, r] += take
                    remaining -= take
                if remaining > 0:
                    raise RuntimeError(
                        f"init_default: Cannot place the population of class {r}: "
                        f"total station capacity is insufficient.")

        # Delegate to initFromMarginal which generates per-node state spaces
        # and state priors (needed for CTMC/FLD transient analysis)
        self.init_from_marginal(n0)

    def has_init_state(self) -> bool:
        """Check if network has initialized state.

        Mirrors MATLAB MNetwork.hasInitState: the model counts as initialized
        only when EVERY stateful node carries a state, so a partially
        initialized model still triggers init_default.
        """
        return all(node.state is not None for node in self._nodes if node.is_stateful())

    def get_state(self) -> list:
        """
        Get current state of all stateful nodes.

        Returns:
            List of state arrays, one per stateful node, in node order.
            Each element is a list containing the node's state array(s).
        """
        if not self.has_init_state():
            self.init_default()
        from .nodes import Fork
        state = []
        for node in self._nodes:
            if node.is_stateful() or (getattr(self, 'fork_stateful', False) and isinstance(node, Fork)):
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

        # store requested marginal flattened [n[0,0],n[0,1],...], padded to nclasses columns, for later CTMC initial-state use.
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

                self._assert_pas_placement_given(sn, station, node_idx, n_row)

                # Generate state space for this station
                # Use fromMarginalAndStarted if started jobs specified, otherwise fromMarginal
                if np.any(s_row > 0):
                    state_space = fromMarginalAndStarted(sn, node_idx, n_row, s_row)
                else:
                    state_space = fromMarginal(sn, node_idx, n_row)

                self._assign_station_state(station, state_space, n_row, s_row)

        self._has_state = True
        self._reset_struct()

    # PascalCase alias
    initFromMarginalAndStarted = init_from_marginal_and_started

    def init_from_marginal_and_running(self, n, s, options=None) -> None:
        """
        Initialize network state from marginal queue lengths and RUNNING jobs.

        Mirrors MATLAB ``@MNetwork/initFromMarginalAndRunning``. It is not
        ``initFromMarginalAndStarted``: 'started' counts the jobs that have
        BEGUN service (so a preempted job still counts), while 'running' counts
        the jobs holding a server right now, which is what a load-dependent or
        multiserver station needs to reproduce a mid-service snapshot.

        Args:
            n: Marginal queue lengths (nstations x nclasses, or nnodes rows);
                n[i][r] is the number of jobs of class r at station i.
            s: Running jobs, same shape as n.
            options: unused, accepted for MATLAB call compatibility.

        Raises:
            ValueError: if the pair (n, s) is not a valid state of this model,
                or if a stateful node ends up with no state.
        """
        from ..api.state.marginal import fromMarginalAndRunning
        from .state import State

        n = np.atleast_2d(np.asarray(n, dtype=float))
        s = np.atleast_2d(np.asarray(s, dtype=float))
        sn = self.getStruct()

        nstations = len(self._stations)
        nnodes = len(self._nodes)
        if nstations < nnodes:
            # A row per STATION or a row per NODE both reach here; the state
            # generator is node-indexed, so a station-indexed matrix is lifted.
            if n.shape[0] == nstations:
                n_nodes = np.zeros((nnodes, n.shape[1]))
                s_nodes = np.zeros((nnodes, s.shape[1]))
                for ist in range(nstations):
                    ind = int(np.asarray(sn.stationToNode).ravel()[ist])
                    n_nodes[ind, :] = n[ist, :]
                    s_nodes[ind, :] = s[ist, :]
                n_by_node, s_by_node = n_nodes, s_nodes
            elif n.shape[0] == nnodes:
                n_by_node, s_by_node = n, s
            else:
                raise ValueError(
                    'The supplied matrix of marginal states does not have the '
                    'correct number of rows. One either per station or per node.')
        else:
            n_by_node, s_by_node = n, s

        if not State.isValid(sn, n, s):
            raise ValueError('Initial state is not valid.')

        nclasses = len(self._classes)
        self._state_marginal = np.asarray(
            [n_by_node[self._nodes.index(st), :nclasses] if st in self._nodes
             else np.zeros(nclasses) for st in self._stations], dtype=float).flatten()

        for i, station in enumerate(self._stations):
            if not hasattr(station, 'set_state'):
                continue
            node_idx = self._nodes.index(station) if station in self._nodes else i
            n_row = n_by_node[node_idx, :nclasses]
            s_row = s_by_node[node_idx, :nclasses]
            self._assert_pas_placement_given(sn, station, node_idx, n_row)
            state_i = fromMarginalAndRunning(sn, node_idx, n_row, s_row)
            if state_i is None or len(np.atleast_2d(state_i)) == 0:
                raise ValueError('Invalid state assignment for station %d' % i)
            self._assign_station_state(station, state_i, n_row, s_row)

        self._has_state = True
        self._reset_struct()

    initFromMarginalAndRunning = init_from_marginal_and_running

    def init_from_avg_qlen(self, AvgQLen) -> None:
        """
        Initialize the state from mean queue lengths, rounded to integers.

        Mirrors MATLAB ``@MNetwork/initFromAvgQLen``. Rounding each entry
        independently can overshoot the closed population -- round([0.5,0.5])
        is [1,1] -- so a class whose rounded total exceeds its mean total gives
        one job back at its fullest station. A marginal that still does not
        validate falls back to the default initialization, as MATLAB's does.

        Args:
            AvgQLen: mean queue lengths (nstations x nclasses).
        """
        AvgQLen = np.atleast_2d(np.asarray(AvgQLen, dtype=float))
        n = np.round(AvgQLen)
        for r in range(AvgQLen.shape[1]):
            # error at most by 1
            if n[:, r].sum() > AvgQLen[:, r].sum():
                i = int(np.argmax(n[:, r]))
                n[i, r] -= 1
        try:
            self.init_from_marginal(n)
        except Exception:
            self.initDefault()

    initFromAvgQLen = init_from_avg_qlen

    def init_from_avg_table_qlen(self, AvgTable) -> None:
        """
        Initialize the state from the QLen column of an average table.

        Mirrors MATLAB ``@MNetwork/initFromAvgTableQLen``: the column is stored
        class-major, so it reshapes to (nclasses, nstations) and transposes.

        Args:
            AvgTable: a DataFrame carrying a 'QLen' column, as getAvgTable returns.
        """
        qlen = np.asarray(AvgTable['QLen'], dtype=float)
        QN = qlen.reshape(len(self._classes), len(self._stations)).T
        self.init_from_avg_qlen(QN)

    initFromAvgTableQLen = init_from_avg_table_qlen

    def init_from_marginal(self, n, options=None) -> None:
        """
        Initialize network state from marginal queue lengths only.

        Mirrors MATLAB ``@MNetwork/initFromMarginal``, which is NOT
        ``initFromMarginalAndStarted`` with a zero started matrix: it validates
        the marginal first, and it admits a PURPOSELY FRACTIONAL one.

        A fractional row is a fluid initial condition, and it is kept as the
        state verbatim. Routing it through the discrete state-space generator
        instead truncates it per station (``astype(int)``), which silently DROPS
        JOBS: SolverENV hands a fluid stage the fractional Qentry, so a closed
        model of 5 jobs was written out holding 4, and MATLAB then refused the
        document with "Chain 1 is initialized with an incorrect number of jobs".

        Args:
            n: Marginal queue lengths. Can be:
               - 1D list/array with one value per station (single class)
               - 2D list/array (nstations x nclasses)
        """
        from ..api.state.marginal import fromMarginal
        from ..constants import GlobalConstants
        from .state import State

        n = np.atleast_2d(np.asarray(n, dtype=float))

        nstations = len(self._stations)
        nclasses = len(self._classes)

        # Handle 1D input (single class - vector of length nstations)
        if n.shape[0] == 1 and nstations > 1:
            # If shape is (1, nstations), transpose to (nstations, 1)
            n = n.T
        if n.shape[1] < nclasses:
            n = np.hstack([n, np.zeros((n.shape[0], nclasses - n.shape[1]))])

        sn = self.getStruct()

        # One recovery attempt, then reject, exactly as MATLAB does: a marginal
        # that is not a state of this network is a caller error, and carrying it
        # silently produces a model whose chain populations are not its own.
        if not State.isValid(sn, n, None):
            n = np.round(n)
            if not State.isValid(sn, n, None):
                raise ValueError(
                    "Initial state not contained in the state space and not "
                    "recoverable by rounding.")

        self._state_marginal = n[:nstations, :nclasses].flatten()

        for i, station in enumerate(self._stations):
            if not hasattr(station, 'set_state'):
                continue
            node_idx = self._nodes.index(station) if station in self._nodes else i
            n_row = n[i, :nclasses]
            self._assert_pas_placement_given(sn, station, node_idx, n_row)
            if np.max(np.abs(n_row - np.round(n_row))) < GlobalConstants.CoarseTol:
                state_space = fromMarginal(sn, node_idx, np.round(n_row))
            else:
                # Purposely fractional, e.g. a fluid solver initialization
                state_space = np.atleast_2d(n_row)
            self._assign_station_state(station, state_space, n_row,
                                       np.zeros(nclasses))

        self._has_state = True
        self._reset_struct()

    def _assert_pas_placement_given(self, sn, station, node_idx, n_row) -> None:
        """A closed PAS station with a nonempty swap graph and >1 initial job has
        a reducible generator; initial placement is a required input, not a
        fabricated default. See _kb/11-conventions-and-gotchas.md (closed PAS
        ordering).
        """
        sched_name = getattr(station.get_sched_strategy(), 'name', None) \
            if hasattr(station, 'get_sched_strategy') else None
        if sched_name != 'PAS' or np.sum(n_row) <= 1:
            return
        sg = None
        if sn.nodeparam is not None and node_idx in sn.nodeparam \
                and isinstance(sn.nodeparam[node_idx], dict):
            sg = sn.nodeparam[node_idx].get('swapGraph')
        has_swap = sg is not None and np.any(np.asarray(sg, dtype=float) != 0)
        is_closed = any(np.isfinite(getattr(jc, '_njobs', np.inf))
                        for jc in self._classes)
        if has_swap and is_closed:
            raise RuntimeError(
                "A closed pass-and-swap station with a non-empty swapping graph "
                "requires an explicit initial job placement. Call setState on the "
                "station with the ordered class list (oldest first) before solving.")

    def _assign_station_state(self, station, state_space, n_row, s_row) -> None:
        """Install a station's state space, its prior and its current state.

        The three are one unit: initDefault, initFromMarginal and
        initFromMarginalAndStarted all set them together, and a node carrying
        any two of them is not a state any solver can start from.
        """
        if state_space is not None and len(state_space) > 0:
            if hasattr(station, 'set_state_space'):
                station.set_state_space(state_space)
            if hasattr(station, 'setStatePrior'):
                prior = np.zeros(len(state_space))
                prior[0] = 1.0
                station.setStatePrior(prior)
            station.set_state(state_space[0])
        else:
            # Fallback to simple state generation
            station.set_state(
                self._state_from_marginal_and_started(station, n_row, s_row))

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
        """Convert Network for JVM interoperability.

        Not available in the native Python implementation. The native Python
        solver operates independently of the JVM. For JVM interoperability call
        the canonical JAR (common/jline.jar) directly.
        """
        raise NotImplementedError(
            "to_java() is not available in the native Python implementation. "
            "The native Python solver operates independently of the JVM. "
            "For JVM interoperability call the canonical JAR (common/jline.jar) directly."
        )

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

        # enableChecks is provenance, not state; see _kb/04-networkstruct.md Python native network.py section.
        new_network._do_checks = getattr(self, '_do_checks', True)

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

            # reference-station lookup after copy must key off the ORIGINAL class's station; see _kb/04-networkstruct.md Python native network.py section.
            if hasattr(old_class, '_refstat') and old_class._refstat is not None:
                if old_class._refstat in node_map:
                    new_class._refstat = node_map[old_class._refstat]
            if hasattr(old_class, '_reference_station') and old_class._reference_station is not None:
                if old_class._reference_station in node_map:
                    new_class._reference_station = node_map[old_class._reference_station]

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

            # deep-copied node: re-key every per-class dict onto the copy's classes; see _kb/04-networkstruct.md Python native network.py section.
            new_node.remap_job_classes(new_network._classes)

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
    getClassNames = get_class_names
    getNodeNames = get_node_names
    getStationNames = get_station_names
    getClassSwitchingMask = get_class_switching_mask
    hasProductFormSolution = has_product_form_solution

    # State management
    refreshStruct = refresh_struct
    resetStruct = reset_struct
    getStruct = get_struct

    # Reward methods
    setReward = set_reward
    getReward = get_reward
    getRewards = get_rewards

    def get_version(self):
        """Get the LINE solver version string."""
        from ..constants import GlobalConstants
        return GlobalConstants.Version

    getVersion = get_version
    version = get_version

    def get_number_of_jobs(self):
        """Total population of the closed classes (finite njobs entries)."""
        total = 0.0
        for c in self._classes:
            n = c.getNumberOfJobs()
            if np.isfinite(n):
                total += n
        return total

    getNumberOfJobs = get_number_of_jobs

    def print_struct(self):
        """Print a summary of the compiled NetworkStruct (wrapper-compatible)."""
        sn = self.get_struct()
        print('nstations: %d  nclasses: %d  nchains: %d  nnodes: %d  nstateful: %d'
              % (sn.nstations, sn.nclasses, sn.nchains, sn.nnodes, sn.nstateful))
        print('nodenames:', list(sn.nodenames))
        print('classnames:', list(sn.classnames))
        print('njobs:', np.asarray(sn.njobs).flatten().tolist())
        print('nservers:', np.asarray(sn.nservers).flatten().tolist())
        print('refstat:', np.asarray(sn.refstat).flatten().tolist())
        schedv = sn.sched
        if isinstance(schedv, dict):
            sched_names = [getattr(schedv[k], 'name', str(schedv[k])) for k in sorted(schedv)]
        else:
            sched_names = [str(x) for x in np.asarray(schedv).flatten().tolist()]
        print('sched:', sched_names)
        print('rates:')
        print(np.asarray(sn.rates))
        print('scv:')
        print(np.asarray(sn.scv))
        if sn.rt is not None:
            print('rt:')
            print(np.asarray(sn.rt))

    printStruct = print_struct

    def get_product_form_parameters(self):
        """Extract product-form parameters (lambda, D, N, Z, mu, S, V), matching
        the wrapper Network.getProductFormParameters return order."""
        from ..api.sn.transforms import sn_get_product_form_params
        params = sn_get_product_form_params(self.get_struct())
        return (params.lam, params.D, params.N, params.Z, params.mu,
                params.S, params.V)

    getProductFormParameters = get_product_form_parameters

    def get_chain_index(self, jobclass):
        """Get the 1-based index of the chain containing a class (object or name)."""
        r = self.get_class_index(jobclass)  # 1-based
        sn = self.get_struct()
        chains = np.asarray(sn.chains)
        for c in range(chains.shape[0]):
            if chains[c, r - 1]:
                return c + 1
        return -1

    getChainIndex = get_chain_index
    chain_index = get_chain_index

    def get_job_class_index(self, jobclass):
        """Get the 1-based index of a job class (alias of get_class_index)."""
        return self.get_class_index(jobclass)

    getJobClassIndex = get_job_class_index

    def get_stateful_index(self, node):
        """Get the 1-based index of a node among the stateful nodes (-1 if not stateful)."""
        sn = self.get_struct()
        ind = node if isinstance(node, int) else self.get_node_index(node)
        mapping = np.asarray(sn.nodeToStateful).flatten()
        if ind - 1 < 0 or ind - 1 >= len(mapping):
            return -1
        isf = int(mapping[ind - 1])
        return isf + 1 if isf >= 0 else -1

    getStatefulIndex = get_stateful_index
    stateful_index = get_stateful_index

    def get_graph(self):
        """Build (H, G) station- and node-level graph dictionaries (see lang/viz.py)."""
        from .viz import network_get_graph
        return network_get_graph(self)

    getGraph = get_graph

    def plot(self, graph_type='station', method='names', **kwargs):
        """Plot the network as a directed graph (see lang/viz.py)."""
        from .viz import network_plot
        return network_plot(self, graph_type=graph_type, method=method, **kwargs)

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
        # serial_routing accepts list, variadic, or numpy-array calling conventions.
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
            Java: jar/src/main/java/jline/lang/Network.java
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
            Java: jar/src/main/java/jline/lang/Network.java
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
            Java: jar/src/main/java/jline/lang/Network.java
        """
        D = np.atleast_2d(D)
        M = D.shape[0]

        # All FCFS scheduling
        strategy = [SchedStrategy.FCFS] * M

        return Network.cyclic(N, D, strategy, S)

    @staticmethod
    def cyclicFcfsInf(N, D, Z=None, S=None):
        """
        Create a cyclic closed network with Delay (INF) stations followed by
        FCFS queues.

        Args:
            N: Population vector [1 x R] or list - number of jobs per class
            D: Service demand matrix [M x R] - demands at FCFS queues per class
            Z: Think time matrix [MZ x R] at delay stations per class
            S: Number of servers per FCFS queue [M x 1] or list (default: 1 each)

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (cyclicFcfsInf)
        """
        D = np.atleast_2d(D)
        M = D.shape[0]
        if Z is None:
            Z = np.zeros((0, D.shape[1]))
        Z = np.atleast_2d(Z)
        MZ = Z.shape[0]
        if MZ > 0 and np.max(Z) == 0:
            MZ = 0
            Z = np.zeros((0, D.shape[1]))

        strategy = [SchedStrategy.INF] * MZ + [SchedStrategy.FCFS] * M
        Dnew = np.vstack([Z, D]) if MZ > 0 else D

        if S is None:
            S = np.ones((M, 1))
        else:
            S = np.atleast_2d(S)
            if S.shape[0] == 1 and S.shape[1] == M:
                S = S.T
        Snew = np.vstack([np.full((MZ, 1), np.inf), S]) if MZ > 0 else S

        return Network.cyclic(N, Dnew, strategy, Snew)

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
            Java: jar/src/main/java/jline/lang/Network.java
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

    def modelView(self) -> bool:
        """
        Open the model in JSIMgraph viewer.

        Exports the network to JSIMG format and launches JMT's JSIMgraph
        as a subprocess to display an interactive visualization.

        Returns:
            True if the viewer was launched successfully, False otherwise

        References:
            MATLAB: matlab/src/lang/@MNetwork/modelView.m
        """
        import tempfile
        import os
        from ..api.io import jsimg_view, line_printf
        # from ..api.io import line_viewer_view
        from ..api.solvers.jmt.handler import _write_jsim_file, SolverJMTOptions

        # Compile model if needed
        if not self._has_struct:
            self.link(self._routing_matrix)

        # Get NetworkStruct
        sn = self.get_struct()

        # Create temp file for JSIMG
        fd, jsimg_file = tempfile.mkstemp(suffix='.jsimg', prefix='model_')
        os.close(fd)

        # Write model to JSIMG format
        options = SolverJMTOptions()
        _write_jsim_file(sn, jsimg_file, options)

        line_printf('JSIMgraph Model: %s\n', jsimg_file)

        # Open in JSIMgraph
        return jsimg_view(jsimg_file)
        # return line_viewer_view(jsimg_file)

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
    def tandemPs(lambda_rates: np.ndarray, D: np.ndarray) -> 'Network':
        """Create an open tandem network of PS queues (no delay stations).

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (tandemPs)
        """
        return Network.tandemPsInf(lambda_rates, D, None)

    @staticmethod
    def tandemFcfsInf(lambda_rates: np.ndarray, D: np.ndarray,
                      Z: Optional[np.ndarray] = None) -> 'Network':
        """Create an open tandem network with FCFS queues and INF (delay)
        stations.

        Args:
            lambda_rates: Array of arrival rates for each class (R,)
            D: Service time matrix for FCFS queues (M x R)
            Z: Optional delay service times (Mz x R)

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (tandemFcfsInf)
        """
        if Z is None:
            Z = np.array([]).reshape(0, D.shape[1]) if len(np.shape(D)) > 1 else np.array([])

        Z = np.atleast_2d(Z) if np.size(Z) > 0 else np.zeros((0, np.atleast_2d(D).shape[1]))
        D = np.atleast_2d(D)

        M = D.shape[0]   # FCFS queues
        Mz = Z.shape[0]  # delay stations

        strategies = [SchedStrategy.INF] * Mz + [SchedStrategy.FCFS] * M
        combined_D = np.vstack([Z, D]) if Mz > 0 else D
        return Network.tandem(lambda_rates, combined_D, strategies)

    @staticmethod
    def tandemFcfs(lambda_rates: np.ndarray, D: np.ndarray) -> 'Network':
        """Create an open tandem network of FCFS queues (no delay stations).

        References:
            MATLAB: matlab/src/lang/@MNetwork/MNetwork.m (tandemFcfs)
        """
        return Network.tandemFcfsInf(lambda_rates, D, None)

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

    @staticmethod
    def cluster(lambda_rates: np.ndarray, D: np.ndarray,
                    strategies: List, S: Optional[np.ndarray] = None,
                    dispatching: 'RoutingStrategy' = None) -> 'Network':
        """
        Create an open server-farm network: Source -> Dispatcher -> Servers -> Sink.

        The dispatcher is a Router that distributes incoming jobs to the M parallel
        server queues according to the supplied dispatching strategy (RAND, RROBIN,
        JSQ, ...).

        Args:
            lambda_rates: Per-class arrival rates (R,)
            D: Service time matrix (M x R); D[i, r] is the mean service time of class
               r at server i
            strategies: Per-server scheduling strategies (length M)
            S: Optional per-server multiplicity (M,) or (M, 1); defaults to all-1
            dispatching: Dispatching policy applied at the router (defaults to RAND)

        Returns:
            Network: Configured open server-farm model
        """
        from .nodes import Source, Queue, Sink, Router
        from .classes import OpenClass
        from ..distributions import Exp

        if dispatching is None:
            dispatching = RoutingStrategy.RAND

        D = np.atleast_2d(D)
        M, R = D.shape

        if S is None:
            S = np.ones(M, dtype=int)
        else:
            S = np.asarray(S).reshape(-1)

        model = Network('Cluster')

        source = Source(model, 'Source')
        dispatcher = Router(model, 'Dispatcher')
        servers = []
        for i in range(M):
            q = Queue(model, f'Station{i+1}', strategies[i])
            if int(S[i]) > 1:
                q.set_number_of_servers(int(S[i]))
            servers.append(q)
        sink = Sink(model, 'Sink')

        lambda_arr = np.atleast_1d(np.asarray(lambda_rates).reshape(-1))
        jobclasses = []
        for r in range(R):
            cls = OpenClass(model, f'Class{r+1}', 0)
            jobclasses.append(cls)
            source.set_arrival(cls, Exp.fit_mean(1.0 / float(lambda_arr[r])))
            for i in range(M):
                servers[i].set_service(cls, Exp.fit_mean(D[i, r]))

        model.add_link(source, dispatcher)
        for q in servers:
            model.add_link(dispatcher, q)
            model.add_link(q, sink)
        for cls in jobclasses:
            dispatcher.set_routing(cls, dispatching)

        return model

    @staticmethod
    def cluster_ps(lambda_rates: np.ndarray, D: np.ndarray,
                       S: Optional[np.ndarray] = None,
                       dispatching: 'RoutingStrategy' = None) -> 'Network':
        """Open PS cluster (single-server queues unless S overrides multiplicity)."""
        D = np.atleast_2d(D)
        M = D.shape[0]
        strategies = [SchedStrategy.PS] * M
        return Network.cluster(lambda_rates, D, strategies, S, dispatching)

    @staticmethod
    def cluster_fcfs(lambda_rates: np.ndarray, D: np.ndarray,
                         S: Optional[np.ndarray] = None,
                         dispatching: 'RoutingStrategy' = None) -> 'Network':
        """Open FCFS cluster."""
        D = np.atleast_2d(D)
        M = D.shape[0]
        strategies = [SchedStrategy.FCFS] * M
        return Network.cluster(lambda_rates, D, strategies, S, dispatching)

    @staticmethod
    def cluster_closed(N: np.ndarray, Z: np.ndarray, D: np.ndarray,
                           strategies: List, S: Optional[np.ndarray] = None,
                           dispatching: 'RoutingStrategy' = None) -> 'Network':
        """
        Create a closed server-farm network: Think -> Dispatcher -> Servers -> Think.

        Args:
            N: Per-class population (1, R) or (R,)
            Z: Per-class think times (1, R) or (R,)
            D: Service time matrix (M, R)
            strategies: Per-server scheduling strategies (length M)
            S: Optional per-server multiplicity (M,)
            dispatching: Dispatching policy (defaults to RAND)
        """
        from .nodes import Queue, Delay, Router
        from .classes import ClosedClass
        from ..distributions import Exp

        if dispatching is None:
            dispatching = RoutingStrategy.RAND

        D = np.atleast_2d(D)
        M, R = D.shape
        N = np.atleast_1d(np.asarray(N).reshape(-1))
        Z = np.atleast_1d(np.asarray(Z).reshape(-1))

        if S is None:
            S = np.ones(M, dtype=int)
        else:
            S = np.asarray(S).reshape(-1)

        model = Network('Cluster')

        think = Delay(model, 'Think')
        dispatcher = Router(model, 'Dispatcher')
        servers = []
        for i in range(M):
            q = Queue(model, f'Station{i+1}', strategies[i])
            if int(S[i]) > 1:
                q.set_number_of_servers(int(S[i]))
            servers.append(q)

        jobclasses = []
        for r in range(R):
            cls = ClosedClass(model, f'Class{r+1}', int(N[r]), think, 0)
            jobclasses.append(cls)
            think.set_service(cls, Exp.fit_mean(float(Z[r])))
            for i in range(M):
                servers[i].set_service(cls, Exp.fit_mean(D[i, r]))

        model.add_link(think, dispatcher)
        for q in servers:
            model.add_link(dispatcher, q)
            model.add_link(q, think)
        for cls in jobclasses:
            dispatcher.set_routing(cls, dispatching)

        return model

    @staticmethod
    def cluster_mixed(lambda_rates: np.ndarray, N: np.ndarray, Z: np.ndarray,
                      D: np.ndarray, strategies: List, S: Optional[np.ndarray] = None,
                      dispatching: 'RoutingStrategy' = None) -> 'Network':
        """
        Create a mixed cluster network sharing dispatcher and servers between
        open and closed classes: open classes flow Source -> Dispatcher ->
        Servers -> Sink, closed classes cycle Think -> Dispatcher -> Servers ->
        Think.

        Classes are ordered open first, so columns 0..Ro-1 of D refer to the
        open classes and columns Ro..Ro+Rc-1 to the closed ones.

        Args:
            lambda_rates: Per-class arrival rates of the open classes (Ro,)
            N: Per-class population of the closed classes (Rc,)
            Z: Per-class think times of the closed classes (Rc,)
            D: Service time matrix (M, Ro+Rc)
            strategies: Per-server scheduling strategies (length M)
            S: Optional per-server multiplicity (M,)
            dispatching: Dispatching policy (defaults to RAND)
        """
        from .nodes import Source, Queue, Sink, Delay, Router
        from .classes import OpenClass, ClosedClass
        from ..distributions import Exp, Disabled

        if dispatching is None:
            dispatching = RoutingStrategy.RAND

        D = np.atleast_2d(D)
        M, R = D.shape
        lambda_arr = np.atleast_1d(np.asarray(lambda_rates, dtype=float).reshape(-1))
        N = np.atleast_1d(np.asarray(N).reshape(-1))
        Z = np.atleast_1d(np.asarray(Z, dtype=float).reshape(-1))
        Ro, Rc = lambda_arr.size, N.size
        if R != Ro + Rc:
            raise ValueError("D must have len(lambda_rates)+len(N) columns")
        if Z.size != Rc:
            raise ValueError("N and Z must have the same length")

        if S is None:
            S = np.ones(M, dtype=int)
        else:
            S = np.asarray(S).reshape(-1)

        model = Network('Cluster')

        source = Source(model, 'Source')
        think = Delay(model, 'Think')
        dispatcher = Router(model, 'Dispatcher')
        servers = []
        for i in range(M):
            q = Queue(model, f'Station{i+1}', strategies[i])
            if int(S[i]) > 1:
                q.set_number_of_servers(int(S[i]))
            servers.append(q)
        sink = Sink(model, 'Sink')

        jobclasses = []
        for r in range(Ro):
            cls = OpenClass(model, f'Class{r+1}', 0)
            jobclasses.append(cls)
            source.set_arrival(cls, Exp.fit_mean(1.0 / float(lambda_arr[r])))
            think.set_service(cls, Disabled())
            for i in range(M):
                servers[i].set_service(cls, Exp.fit_mean(D[i, r]))
        for c in range(Rc):
            r = Ro + c
            cls = ClosedClass(model, f'Class{r+1}', int(N[c]), think, 0)
            jobclasses.append(cls)
            think.set_service(cls, Exp.fit_mean(float(Z[c])))
            for i in range(M):
                servers[i].set_service(cls, Exp.fit_mean(D[i, r]))

        model.add_link(source, dispatcher)
        model.add_link(think, dispatcher)
        for q in servers:
            model.add_link(dispatcher, q)
            model.add_link(q, sink)
            model.add_link(q, think)

        # Class-specific exits: the servers feed the sink for open classes and
        # the delay for closed ones, so the shared arcs carry explicit
        # per-class probabilities.
        for r, cls in enumerate(jobclasses):
            dispatcher.set_routing(cls, dispatching)
            is_open = 1.0 if r < Ro else 0.0
            for q in servers:
                q.set_prob_routing(cls, sink, is_open)
                q.set_prob_routing(cls, think, 1.0 - is_open)
            source.set_prob_routing(cls, dispatcher, is_open)
            think.set_prob_routing(cls, dispatcher, 1.0 - is_open)

        return model

    # Snake_case aliases
    tandem_ps_inf = tandemPsInf
    tandem_ps = tandemPs
    tandem_fcfs = tandemFcfs
    tandem_fcfs_inf = tandemFcfsInf
    cyclic_ps_inf = cyclicPsInf
    cyclic_fcfs = cyclicFcfs
    cyclic_fcfs_inf = cyclicFcfsInf
    cyclic_ps = cyclicPs

    # CamelCase aliases for the server-farm factories (mirror the Java naming)
    cluster = cluster
    clusterPs = cluster_ps
    clusterFcfs = cluster_fcfs
    clusterClosed = cluster_closed
    clusterMixed = cluster_mixed


__all__ = ['Network']


def _sdr_same(a, b):
    """Structural equality of two SDR declarations, ignoring the entry node."""
    import numpy as _np
    keys = ('departure', 'branch', 'entryOf', 'departureOf', 'level', 'C')
    for k in keys:
        if a.get(k) != b.get(k):
            return False
    return _np.array_equal(_np.asarray(a['d']), _np.asarray(b['d']))
