"""
Native Python implementation of NetworkGenerator.

Generates random queueing network models without requiring the Java backend.
This is a pure Python implementation of the random network generation functionality.
"""

import random
import numpy as np
from typing import Callable, List, Optional, Tuple

from ..lang.base import RoutingStrategy, SchedStrategy
from ..lang.classes import ClosedClass, OpenClass
from ..lang.network import Network
from ..lang.nodes import ClassSwitch, Delay, Queue, Sink, Source
from ..distributions import Erlang, Exp, HyperExp


# Constants
MAX_SERVERS = 40
HIGH_JOB_LOAD_RANGE = (31, 40)
MED_JOB_LOAD_RANGE = (11, 20)
LOW_JOB_LOAD_RANGE = (1, 5)


def rand_spanning_tree(num_vertices: int) -> np.ndarray:
    """
    Generate a random spanning tree with the specified number of vertices.

    Creates a random tree structure by connecting each vertex i (i > 0)
    to a random vertex j where j < i.

    Args:
        num_vertices: Number of vertices in the tree.

    Returns:
        Adjacency matrix of shape (num_vertices, num_vertices).
    """
    tree = np.zeros((num_vertices, num_vertices))
    for i in range(1, num_vertices):
        parent = random.randint(0, i - 1)
        tree[parent, i] = 1.0
    return tree


def rand_graph(num_vertices: int) -> np.ndarray:
    """
    Generate a random strongly connected graph topology.

    Creates a random strongly connected directed graph with the specified
    number of vertices using a modified Tarjan's algorithm approach.
    The algorithm ensures all vertices are reachable from each other.

    Args:
        num_vertices: Number of vertices in the graph.

    Returns:
        Adjacency matrix of shape (num_vertices, num_vertices).

    Raises:
        ValueError: If num_vertices is not positive.

    Examples:
        >>> adj = rand_graph(4)
        >>> adj.shape
        (4, 4)
    """
    if num_vertices <= 0:
        raise ValueError("Number of vertices must be positive")

    if num_vertices == 1:
        adj = np.zeros((1, 1))
        adj[0, 0] = 1.0
        return adj

    # State variables for DFS-based strong connectivity
    global_start_time = [0]  # Use list for mutable closure
    start_time = np.full(num_vertices, -1)
    lowest_link = np.full(num_vertices, -1)
    inv_start_time = np.full(num_vertices, -1)

    # Generate random spanning tree
    tree = rand_spanning_tree(num_vertices)

    def strong_connect(g: np.ndarray, v: int) -> np.ndarray:
        """DFS-based algorithm to ensure strong connectivity."""
        start_time[v] = global_start_time[0]
        inv_start_time[start_time[v]] = v
        lowest_link[v] = start_time[v]
        global_start_time[0] += 1

        # Find outgoing edges from v
        out_vertices = np.where(g[v, :] > 0)[0]

        for w in out_vertices:
            if start_time[w] == -1:
                g = strong_connect(g, w)
                lowest_link[v] = min(lowest_link[v], lowest_link[w])
            else:
                lowest_link[v] = min(lowest_link[v], start_time[w])

        # If v is root of SCC but not the entire graph's root
        if lowest_link[v] == start_time[v] and start_time[v] > 0:
            descendant_st = random.randint(start_time[v], global_start_time[0] - 1)
            ancestor_st = random.randint(0, start_time[v] - 1)
            # Add edge to ensure strong connectivity
            g[inv_start_time[descendant_st], inv_start_time[ancestor_st]] = 1.0
            lowest_link[v] = ancestor_st

        return g

    # Apply strong connectivity algorithm
    strong_graph = strong_connect(tree.copy(), 0)

    # Randomly permute nodes
    perm = list(range(num_vertices))
    random.shuffle(perm)
    permuted = np.zeros((num_vertices, num_vertices))
    for i in range(num_vertices):
        for j in range(num_vertices):
            if strong_graph[i, j] > 0:
                permuted[perm[i], perm[j]] = 1.0

    return permuted


def cyclic_graph(num_vertices: int) -> np.ndarray:
    """
    Generate a cyclic graph topology.

    Creates a cyclic directed graph where each vertex i is connected
    to vertex (i+1) mod n.

    Args:
        num_vertices: Number of vertices in the graph.

    Returns:
        Adjacency matrix of shape (num_vertices, num_vertices).

    Raises:
        ValueError: If num_vertices is not positive.

    Examples:
        >>> adj = cyclic_graph(4)
        >>> # Creates edges: 0->1, 1->2, 2->3, 3->0
    """
    if num_vertices <= 0:
        raise ValueError("Number of vertices must be positive")

    adj = np.zeros((num_vertices, num_vertices))

    if num_vertices == 1:
        adj[0, 0] = 1.0
    else:
        for i in range(num_vertices - 1):
            adj[i, i + 1] = 1.0
        adj[num_vertices - 1, 0] = 1.0

    return adj


def randfixedsumone(num_elems: int) -> List[float]:
    """
    Generate random probabilities that sum to 1.

    Args:
        num_elems: Number of probability elements.

    Returns:
        List of probabilities summing to 1.
    """
    if num_elems == 0:
        return []
    if num_elems == 1:
        return [1.0]

    # Generate random values and normalize
    values = [random.random() for _ in range(num_elems)]
    total = sum(values)

    # Normalize to sum to 1
    probs = [np.ceil(v / total * 1000) / 1000 for v in values]

    # Adjust largest element to ensure exact sum of 1
    max_idx = probs.index(max(probs))
    current_sum = sum(probs)
    probs[max_idx] -= (current_sum - 1.0)

    return probs


def randintfixedsum(s: int, n: int) -> List[int]:
    """
    Generate n random positive integers that sum to s.

    Args:
        s: Target sum.
        n: Number of integers.

    Returns:
        List of positive integers summing to s.
    """
    if n == 1:
        return [s]
    elif s == n:
        return [1] * n

    first = random.randint(1, s - n)
    rest = randintfixedsum(s - first, n - 1)
    result = [first] + rest
    random.shuffle(result)
    return result


class NetworkGenerator:
    """
    A pure Python generator for creating random queueing network models.

    This class generates random queueing networks according to the
    configurable properties below. It does not depend on the Java backend.

    Two entry points are provided:
        generate: returns a fully constructed Network object. This mirrors
            the MATLAB NetworkGenerator.generate and the Java
            jline.gen.NetworkGenerator.generate.
        generate_config: returns the same random draws as a plain dictionary,
            for callers that want the configuration without a Network object.

    Attributes:
        sched_strat: Scheduling strategy ('fcfs', 'ps', or 'randomize').
        routing_strat: Routing strategy ('Probabilities', 'Random', or 'randomize').
        distribution: Service distribution ('Exp', 'Erlang', 'HyperExp', or 'randomize').
        cclass_job_load: Closed class job load ('low', 'medium', 'high', or 'randomize').
        has_varying_service_rates: Whether to vary service rates.
        has_multi_server_queues: Whether to allow multi-server queues.
        has_random_cs_nodes: Whether to add random class switch nodes.
        has_multi_chain_cs: Whether to allow multi-chain class switching.
        topology_fcn: Function to generate network topology.

    Examples:
        >>> gen = NetworkGenerator()
        >>> model = gen.generate(3, 1, 1, 2)
        >>> config = gen.generate_config(3, 1, 1, 2)
    """

    VALID_SCHED_STRATS = {'fcfs', 'ps', 'inf', 'lcfs', 'lcfspr', 'siro',
                          'sjf', 'ljf', 'sept', 'lept', 'randomize'}
    VALID_ROUTING_STRATS = {'Probabilities', 'Random', 'randomize'}
    VALID_DISTRIBUTIONS = {'exp', 'erlang', 'hyperexp', 'randomize'}
    VALID_JOB_LOADS = {'low', 'medium', 'high', 'randomize'}

    def __init__(
        self,
        sched_strat: str = 'randomize',
        routing_strat: str = 'randomize',
        distribution: str = 'randomize',
        cclass_job_load: str = 'randomize',
        has_varying_service_rates: bool = True,
        has_multi_server_queues: bool = True,
        has_random_cs_nodes: bool = True,
        has_multi_chain_cs: bool = True,
        topology_fcn: Optional[Callable[[int], np.ndarray]] = None
    ):
        """
        Initialize a NetworkGenerator with configurable properties.

        Args:
            sched_strat: Scheduling strategy. Options: 'fcfs', 'ps', 'inf',
                'lcfs', 'lcfspr', 'siro', 'sjf', 'ljf', 'sept', 'lept', 'randomize'.
            routing_strat: Routing strategy. Options: 'Probabilities',
                'Random', 'randomize'.
            distribution: Service distribution. Options: 'Exp', 'Erlang',
                'HyperExp', 'randomize'.
            cclass_job_load: Closed class job load. Options: 'low',
                'medium', 'high', 'randomize'.
            has_varying_service_rates: Whether to vary service rates.
            has_multi_server_queues: Whether to allow multi-server queues.
            has_random_cs_nodes: Whether to add random class switch nodes.
            has_multi_chain_cs: Whether to allow multi-chain class switching.
            topology_fcn: Function that takes an integer and returns an
                adjacency matrix. Default: rand_graph.
        """
        self.sched_strat = sched_strat
        self.routing_strat = routing_strat
        self.distribution = distribution
        self.cclass_job_load = cclass_job_load
        self.has_varying_service_rates = has_varying_service_rates
        self.has_multi_server_queues = has_multi_server_queues
        self.has_random_cs_nodes = has_random_cs_nodes
        self.has_multi_chain_cs = has_multi_chain_cs
        self.topology_fcn = topology_fcn if topology_fcn is not None else rand_graph

    @property
    def sched_strat(self) -> str:
        """Get the scheduling strategy."""
        return self._sched_strat

    @sched_strat.setter
    def sched_strat(self, value: str):
        """Set the scheduling strategy with validation."""
        if value.lower() not in self.VALID_SCHED_STRATS:
            raise ValueError(f"Scheduling strategy '{value}' not supported. "
                           f"Valid options: {self.VALID_SCHED_STRATS}")
        self._sched_strat = value.lower()

    @property
    def routing_strat(self) -> str:
        """Get the routing strategy."""
        return self._routing_strat

    @routing_strat.setter
    def routing_strat(self, value: str):
        """Set the routing strategy with validation."""
        if value not in self.VALID_ROUTING_STRATS:
            raise ValueError(f"Routing strategy '{value}' not supported. "
                           f"Valid options: {self.VALID_ROUTING_STRATS}")
        self._routing_strat = value

    @property
    def distribution(self) -> str:
        """Get the distribution type."""
        return self._distribution

    @distribution.setter
    def distribution(self, value: str):
        """Set the distribution type with validation."""
        if value.lower() not in self.VALID_DISTRIBUTIONS:
            raise ValueError(f"Distribution '{value}' not supported. "
                           f"Valid options: {self.VALID_DISTRIBUTIONS}")
        self._distribution = value.lower()

    @property
    def cclass_job_load(self) -> str:
        """Get the closed class job load setting."""
        return self._cclass_job_load

    @cclass_job_load.setter
    def cclass_job_load(self, value: str):
        """Set the closed class job load with validation."""
        if value.lower() not in self.VALID_JOB_LOADS:
            raise ValueError(f"Job load '{value}' not supported. "
                           f"Valid options: {self.VALID_JOB_LOADS}")
        self._cclass_job_load = value.lower()

    def _choose_sched_strat(self) -> str:
        """Choose a scheduling strategy based on settings."""
        if self.sched_strat == 'randomize':
            return random.choice(['fcfs', 'ps'])
        return self.sched_strat

    def _choose_routing_strat(self) -> str:
        """Choose a routing strategy based on settings."""
        if self.routing_strat == 'randomize':
            return random.choice(['Random', 'Probabilities'])
        return self.routing_strat

    def _choose_distribution(self) -> Tuple[str, dict]:
        """
        Choose a distribution and its parameters.

        Returns:
            Tuple of (distribution_name, parameters_dict).
        """
        if self.distribution == 'randomize':
            dist_type = random.choice(['exp', 'erlang', 'hyperexp'])
        else:
            dist_type = self.distribution

        # Choose mean
        if self.has_varying_service_rates:
            mean = 2.0 ** random.randint(-6, 6)
        else:
            mean = 1.0

        if dist_type == 'exp':
            return ('exp', {'rate': 1.0 / mean})
        elif dist_type == 'erlang':
            k = 2 ** random.randint(0, 6)
            return ('erlang', {'rate': k / mean, 'shape': k})
        else:  # hyperexp
            scv = 2.0 ** random.randint(0, 6)
            return ('hyperexp', {'mean': mean, 'scv': scv})

    def _choose_num_jobs(self) -> int:
        """Choose number of jobs based on load settings."""
        if self.cclass_job_load == 'high':
            return random.randint(*HIGH_JOB_LOAD_RANGE)
        elif self.cclass_job_load == 'medium':
            return random.randint(*MED_JOB_LOAD_RANGE)
        elif self.cclass_job_load == 'low':
            return random.randint(*LOW_JOB_LOAD_RANGE)
        else:  # randomize
            return random.randint(1, HIGH_JOB_LOAD_RANGE[1])

    def _choose_num_servers(self) -> int:
        """Choose number of servers based on settings."""
        if self.has_multi_server_queues:
            return random.randint(1, MAX_SERVERS)
        return 1

    def _generate_cs_mask(self, num_oclass: int, num_cclass: int) -> np.ndarray:
        """
        Generate class switching mask.

        Args:
            num_oclass: Number of open classes.
            num_cclass: Number of closed classes.

        Returns:
            Boolean mask matrix indicating valid class switching pairs.
        """
        total_classes = num_oclass + num_cclass
        mask = np.zeros((total_classes, total_classes))

        if not self.has_multi_chain_cs:
            # Open classes can switch among themselves
            mask[:num_oclass, :num_oclass] = 1.0
            # Closed classes can switch among themselves
            mask[num_oclass:, num_oclass:] = 1.0
            return mask

        # Multi-chain class switching logic
        all_chains = []
        if num_oclass > 0:
            num_open_chains = random.randint(1, num_oclass)
            all_chains.extend(randintfixedsum(num_oclass, num_open_chains))
        if num_cclass > 0:
            num_closed_chains = random.randint(1, num_cclass)
            all_chains.extend(randintfixedsum(num_cclass, num_closed_chains))

        start_idx = 0
        for chain_size in all_chains:
            end_idx = start_idx + chain_size
            mask[start_idx:end_idx, start_idx:end_idx] = 1.0
            start_idx = end_idx

        return mask

    def generate_config(
        self,
        num_queues: Optional[int] = None,
        num_delays: Optional[int] = None,
        num_oclass: int = 0,
        num_cclass: Optional[int] = None
    ) -> dict:
        """
        Generate a random network configuration.

        Args:
            num_queues: Number of queues. If None, randomly chosen (1-8).
            num_delays: Number of delay nodes. If None, randomly chosen.
            num_oclass: Number of open classes. Default: 0.
            num_cclass: Number of closed classes. If None, randomly chosen (1-4).

        Returns:
            Dictionary containing network configuration:
                - num_queues: Number of queues
                - num_delays: Number of delays
                - num_oclass: Number of open classes
                - num_cclass: Number of closed classes
                - queues: List of queue configurations
                - delays: List of delay configurations
                - open_classes: List of open class configurations
                - closed_classes: List of closed class configurations
                - service_processes: Service distribution for each (station, class) pair
                - topology: Adjacency matrix for station connectivity
                - routing: Routing configuration per station/class
                - class_switches: List of class switch node configurations
                - cs_mask: Class switching mask matrix

        Raises:
            ValueError: If arguments are invalid.

        Examples:
            >>> gen = NetworkGenerator()
            >>> config = gen.generate_config(3, 1, 1, 2)
        """
        # Set defaults
        if num_queues is None:
            num_queues = random.randint(1, 8)
        if num_delays is None:
            if num_queues > 1:
                num_delays = random.randint(0, 1)
            else:
                num_delays = 1
        if num_cclass is None:
            num_cclass = random.randint(1, 4)

        # Validate
        if num_queues < 0 or num_delays < 0 or num_oclass < 0 or num_cclass < 0:
            raise ValueError("Arguments must be non-negative")
        if num_queues + num_delays <= 0:
            raise ValueError("At least one station required")
        if num_oclass + num_cclass <= 0:
            raise ValueError("At least one job class required")

        num_stations = num_queues + num_delays

        # Create queue configurations
        queues = []
        for i in range(num_queues):
            queues.append({
                'name': f'queue{i + 1}',
                'sched_strat': self._choose_sched_strat(),
                'num_servers': self._choose_num_servers()
            })

        # Create delay configurations
        delays = []
        for i in range(num_delays):
            delays.append({
                'name': f'delay{i + 1}'
            })

        # Create open class configurations
        open_classes = []
        for i in range(num_oclass):
            dist_type, dist_params = self._choose_distribution()
            open_classes.append({
                'name': f'OClass{i + 1}',
                'arrival_dist': dist_type,
                'arrival_params': dist_params
            })

        # Create closed class configurations
        closed_classes = []
        ref_station_idx = random.randint(0, num_stations - 1)
        for i in range(num_cclass):
            closed_classes.append({
                'name': f'CClass{i + 1}',
                'num_jobs': self._choose_num_jobs(),
                'ref_station_idx': ref_station_idx
            })

        # Create service processes
        service_processes = {}
        total_classes = num_oclass + num_cclass
        for station_idx in range(num_stations):
            for class_idx in range(total_classes):
                dist_type, dist_params = self._choose_distribution()
                service_processes[(station_idx, class_idx)] = {
                    'dist_type': dist_type,
                    'params': dist_params
                }

        # Generate topology
        topology = self.topology_fcn(num_stations)

        # Generate class switching mask
        cs_mask = self._generate_cs_mask(num_oclass, num_cclass)

        # Extend topology for source/sink if open classes exist
        if num_oclass > 0:
            # Add source (index num_stations) and sink (index num_stations + 1)
            extended_topology = np.zeros((num_stations + 2, num_stations + 2))
            extended_topology[:num_stations, :num_stations] = topology
            # Source connects to random station
            source_dest = random.randint(0, num_stations - 1)
            extended_topology[num_stations, source_dest] = 1.0
            # Random station connects to sink
            sink_source = random.randint(0, num_stations - 1)
            extended_topology[sink_source, num_stations + 1] = 1.0
            topology = extended_topology

        # Generate routing configurations
        routing = {}
        class_switches = []
        for station_idx in range(num_stations):
            dest_indices = np.where(topology[station_idx, :num_stations] > 0)[0].tolist()
            station_routing = {}

            for class_idx in range(total_classes):
                routing_strat = self._choose_routing_strat()
                station_routing[class_idx] = {
                    'strategy': routing_strat,
                    'destinations': dest_indices
                }

                if routing_strat == 'Probabilities' and dest_indices:
                    # Generate routing probabilities
                    probs = randfixedsumone(len(dest_indices))
                    station_routing[class_idx]['probabilities'] = dict(zip(dest_indices, probs))

            routing[station_idx] = station_routing

            # Potentially add class switch nodes
            if self.has_random_cs_nodes:
                for dest_idx in dest_indices:
                    if random.choice([True, False]):
                        # Generate class switch matrix
                        cs_matrix = np.zeros((total_classes, total_classes))
                        for i in range(total_classes):
                            valid_targets = np.where(cs_mask[i, :] > 0)[0]
                            if len(valid_targets) > 0:
                                probs = randfixedsumone(len(valid_targets))
                                for k, j in enumerate(valid_targets):
                                    cs_matrix[i, j] = probs[k]

                        class_switches.append({
                            'name': f'cs_{station_idx}_{dest_idx}',
                            'source_station': station_idx,
                            'dest_station': dest_idx,
                            'cs_matrix': cs_matrix
                        })

        return {
            'num_queues': num_queues,
            'num_delays': num_delays,
            'num_oclass': num_oclass,
            'num_cclass': num_cclass,
            'queues': queues,
            'delays': delays,
            'open_classes': open_classes,
            'closed_classes': closed_classes,
            'service_processes': service_processes,
            'topology': topology,
            'routing': routing,
            'class_switches': class_switches,
            'cs_mask': cs_mask
        }

    def _make_distribution(self):
        """Draw a random service or arrival distribution object."""
        dist_type, params = self._choose_distribution()
        if dist_type == 'exp':
            return Exp(params['rate'])
        elif dist_type == 'erlang':
            return Erlang(params['rate'], params['shape'])
        else:
            return HyperExp.fitMeanAndSCV(params['mean'], params['scv'])

    @staticmethod
    def _validate_args(num_queues, num_delays, num_oclass, num_cclass):
        """Validate that parameter values for the network are sound."""
        if num_queues < 0 or num_delays < 0 or num_oclass < 0 or num_cclass < 0:
            raise ValueError("Arguments must be non-negative")
        if num_queues + num_delays <= 0 or num_oclass + num_cclass <= 0:
            raise ValueError("At least one station and one job class required")

    def _create_stations(self, model, num_queues, num_delays, has_oclass):
        """Create the stations in the network, including source and sink.

        Source and sink are intentionally excluded from the returned list:
        service processes are not set on them by the caller.
        """
        stations = []
        for i in range(num_queues):
            queue = Queue(model, f'queue{i + 1}', self._choose_sched_strat_enum())
            queue.set_number_of_servers(self._choose_num_servers())
            stations.append(queue)

        for i in range(num_delays):
            stations.append(Delay(model, f'delay{i + 1}'))

        if has_oclass:
            Source(model, 'source')
            Sink(model, 'sink')

        return stations

    def _create_classes(self, model, num_oclass, num_cclass):
        """Create the job classes serviced within the network.

        Arrival processes for open classes are set here, mirroring MATLAB.
        """
        classes = []
        for i in range(num_oclass):
            open_class = OpenClass(model, f'OClass{i + 1}')
            model.get_source().set_arrival(open_class, self._make_distribution())
            classes.append(open_class)

        # The reference station is drawn among the non-Source stations
        num_refstat_candidates = model.get_number_of_stations() - (1 if num_oclass > 0 else 0)
        refstat = model.get_stations()[random.randint(0, num_refstat_candidates - 1)]
        for i in range(num_cclass):
            classes.append(ClosedClass(model, f'CClass{i + 1}',
                                       self._choose_num_jobs(), refstat))

        return classes

    def _set_service_processes(self, stations, classes):
        """Set the service processes for all job classes at all stations."""
        for jobclass in classes:
            for station in stations:
                station.set_service(jobclass, self._make_distribution())

    def _define_topology(self, model, topology):
        """Define the topology of the network, adding random class switch nodes."""
        num_stations = model.get_number_of_stations()
        cs_mask = self._gen_cs_mask(model)

        if model.has_open_classes():
            # Extend the station-indexed adjacency matrix with source and sink
            num_nodes = model.get_number_of_nodes()
            extended = np.zeros((num_nodes, num_nodes))
            extended[:topology.shape[0], :topology.shape[1]] = topology
            extended[model.get_index_source_node(), random.randint(0, num_stations - 2)] = 1.0
            extended[random.randint(0, num_stations - 2), model.get_index_sink_node()] = 1.0
            topology = extended

        for i in range(num_stations):
            dest_ids = sorted(np.where(topology[i, :] > 0)[0].tolist())
            outgoing_nodes = self._add_outgoing_links(model, i, dest_ids, cs_mask)
            self._set_routing_strategies(model, model.get_stations()[i], outgoing_nodes)

    def _gen_cs_mask(self, model):
        """Create a mask indicating which classes can switch with each other."""
        num_oclass = sum(1 for c in model.get_classes() if isinstance(c, OpenClass))
        num_cclass = sum(1 for c in model.get_classes() if isinstance(c, ClosedClass))
        return self._generate_cs_mask(num_oclass, num_cclass)

    def _add_outgoing_links(self, model, source_id, all_dest_ids, mask):
        """For a single node, add all outgoing links and random class switch nodes."""
        source_node = model.get_nodes()[source_id]
        outgoing_nodes = []

        for dest_id in all_dest_ids:
            dest_node = model.get_nodes()[dest_id]
            if (self.has_random_cs_nodes
                    and random.randint(0, 1) == 1
                    and not isinstance(source_node, Source)
                    and not isinstance(dest_node, Sink)):
                cs = self._rand_class_switch_node(model, mask, source_node, dest_node)
                outgoing_nodes.append(cs)
                model.add_link(source_node, cs)
                model.add_link(cs, dest_node)
                for jobclass in model.get_classes():
                    cs.set_prob_routing(jobclass, dest_node, 1.0)
            else:
                outgoing_nodes.append(dest_node)
                model.add_link(source_node, dest_node)

        return outgoing_nodes

    def _rand_class_switch_node(self, model, mask, source_node, dest_node):
        """Instantiate a class switch node with a random switching matrix."""
        name = f'cs_{source_node.get_name()}_{dest_node.get_name()}'
        size = mask.shape[0]
        matrix = np.zeros((size, size))
        for i in range(size):
            targets = np.where(mask[i, :] > 0)[0]
            if len(targets) > 0:
                probs = randfixedsumone(len(targets))
                for k, j in enumerate(targets):
                    matrix[i, j] = probs[k]
        return ClassSwitch(model, name, matrix)

    def _set_routing_strategies(self, model, station, outgoing_nodes):
        """Set routing strategies for each job class at a specified station."""
        if isinstance(station, Source):
            for jobclass in model.get_classes():
                if isinstance(jobclass, OpenClass):
                    for node in outgoing_nodes:
                        station.set_prob_routing(jobclass, node, 1.0)
            return

        for jobclass in model.get_classes():
            strat = self._choose_routing_strat()
            station.set_routing(jobclass, RoutingStrategy.PROB
                                if strat == 'Probabilities' else RoutingStrategy.RAND)
            if strat == 'Probabilities':
                # Closed classes cannot route to a sink
                if isinstance(jobclass, ClosedClass) and outgoing_nodes and \
                        isinstance(outgoing_nodes[-1], Sink):
                    station.set_prob_routing(jobclass, outgoing_nodes[-1], 0.0)
                    probs = randfixedsumone(len(outgoing_nodes) - 1)
                else:
                    probs = randfixedsumone(len(outgoing_nodes))
                for j, prob in enumerate(probs):
                    station.set_prob_routing(jobclass, outgoing_nodes[j], prob)

    def _choose_sched_strat_enum(self) -> SchedStrategy:
        """Choose a scheduling strategy as a SchedStrategy enum member."""
        return getattr(SchedStrategy, self._choose_sched_strat().upper())

    def generate(
        self,
        num_queues: Optional[int] = None,
        num_delays: Optional[int] = None,
        num_oclass: int = 0,
        num_cclass: Optional[int] = None
    ) -> Network:
        """
        Generate a random queueing network model.

        Returns a Network built according to the generator's properties,
        mirroring MATLAB's NetworkGenerator.generate.

        Args:
            num_queues: Number of queues. If None, randomly chosen (1-8).
            num_delays: Number of delay nodes. If None, randomly chosen.
            num_oclass: Number of open classes. Default: 0.
            num_cclass: Number of closed classes. If None, randomly chosen (1-4).

        Returns:
            A fully constructed Network object.

        Raises:
            ValueError: If arguments are invalid.

        Examples:
            >>> gen = NetworkGenerator()
            >>> model = gen.generate(3, 1, 1, 2)
        """
        if num_queues is None:
            num_queues = random.randint(1, 8)
        if num_delays is None:
            num_delays = random.randint(0, 1) if num_queues > 1 else 1
        if num_cclass is None:
            num_cclass = random.randint(1, 4)

        self._validate_args(num_queues, num_delays, num_oclass, num_cclass)

        model = Network('nw')
        stations = self._create_stations(model, num_queues, num_delays, num_oclass > 0)
        classes = self._create_classes(model, num_oclass, num_cclass)
        self._set_service_processes(stations, classes)
        self._define_topology(model, self.topology_fcn(len(stations)))
        return model

    # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
    @property
    def schedStrat(self): return self.sched_strat
    @schedStrat.setter
    def schedStrat(self, v): self.sched_strat = v

    @property
    def routingStrat(self): return self.routing_strat
    @routingStrat.setter
    def routingStrat(self, v): self.routing_strat = v

    @property
    def cclassJobLoad(self): return self.cclass_job_load
    @cclassJobLoad.setter
    def cclassJobLoad(self, v): self.cclass_job_load = v

    @property
    def hasVaryingServiceRates(self): return self.has_varying_service_rates
    @hasVaryingServiceRates.setter
    def hasVaryingServiceRates(self, v): self.has_varying_service_rates = v

    @property
    def hasMultiServerQueues(self): return self.has_multi_server_queues
    @hasMultiServerQueues.setter
    def hasMultiServerQueues(self, v): self.has_multi_server_queues = v

    @property
    def hasRandomCSNodes(self): return self.has_random_cs_nodes
    @hasRandomCSNodes.setter
    def hasRandomCSNodes(self, v): self.has_random_cs_nodes = v

    @property
    def hasMultiChainCS(self): return self.has_multi_chain_cs
    @hasMultiChainCS.setter
    def hasMultiChainCS(self, v): self.has_multi_chain_cs = v

    randGraph = staticmethod(rand_graph)
    cyclicGraph = staticmethod(cyclic_graph)
