"""
JMT Solver handler - Native Python implementation.

Calls JMT via subprocess.

Port from:



Note: This is a simplified implementation supporting basic queueing networks.
Complex features (caches, transitions, etc.) may require the Java implementation.
"""

import numpy as np
import subprocess
import tempfile
import os
import shutil
import platform
import urllib.request
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
from xml.etree import ElementTree as ET
from xml.dom import minidom
import time

from ...sn import (
    NetworkStruct,
    NodeType,
    SchedStrategy,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
)
from ....constants import ProcessType


@dataclass
class SolverJMTOptions:
    """Options for JMT solver."""
    method: str = 'jsim'
    samples: int = 10000
    seed: int = 23000
    max_simulated_time: float = float('inf')
    conf_int: float = 0.99
    max_rel_err: float = 0.03
    verbose: bool = False


@dataclass
class SolverJMTReturn:
    """
    Result of JMT solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        A: Arrival rates (M x K)
        W: Waiting times (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        runtime: Runtime in seconds
        method: Method used
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    W: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "jsim"


def _get_jmt_jar_path() -> str:
    """Get path to JMT.jar, downloading if necessary."""
    # Look in common/ directory
    package_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    python_dir = os.path.dirname(package_dir)
    root_dir = os.path.dirname(python_dir)
    common_dir = os.path.join(root_dir, 'common')
    jmt_path = os.path.join(common_dir, 'JMT.jar')

    if os.path.isfile(jmt_path):
        return jmt_path

    # Try to download
    os.makedirs(common_dir, exist_ok=True)
    jmt_url = 'https://line-solver.sourceforge.net/latest/JMT.jar'
    try:
        urllib.request.urlretrieve(jmt_url, jmt_path)
        return jmt_path
    except Exception as e:
        raise RuntimeError(
            f"JMT.jar not found and download failed: {e}\n"
            f"Please manually download from {jmt_url} and place in {common_dir}"
        )


def is_jmt_available() -> bool:
    """Check if JMT is available."""
    # Check for JMT
    try:
        result = subprocess.run(
            ['java', '-version'],
            capture_output=True,
            timeout=5
        )
        if result.returncode != 0:
            return False
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False

    # Check for JMT.jar
    try:
        jmt_path = _get_jmt_jar_path()
        return os.path.isfile(jmt_path)
    except RuntimeError:
        return False


def _get_sched_strategy_class(sched: SchedStrategy) -> str:
    """Map scheduling strategy to JMT QueueGetStrategy class name.

    Note: JMT handles PS/GPS/DPS through PSStrategies, not QueueGetStrategies.
    For QueueGetStrategies, only FCFS and LCFS are available.
    PS and other strategies are handled by the Server section's PSStrategy.
    """
    strategy_map = {
        SchedStrategy.FCFS: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.LCFS: "jmt.engine.NetStrategies.QueueGetStrategies.LCFSstrategy",
        # PS/SIRO/INF use FCFS queue get strategy - actual scheduling is in PSStrategy
        SchedStrategy.PS: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.SIRO: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.INF: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
    }
    return strategy_map.get(sched, "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy")


def _detect_class_switches(sn: NetworkStruct) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Detect node pairs that require ClassSwitch nodes for JMT.

    For each pair of connected nodes (i, j), check if there's class switching
    in the routing. If the class transition matrix from i to j is not diagonal,
    a ClassSwitch node needs to be inserted.

    Args:
        sn: NetworkStruct object

    Returns:
        Dictionary mapping (source_idx, dest_idx) to class switching matrix (K x K)
        where matrix[r][s] is the probability of switching from class r to class s
        when routing from source to dest.
    """
    if sn.rtnodes is None or sn.connmatrix is None:
        return {}

    K = sn.nclasses
    nnodes = sn.nnodes
    cs_nodes = {}

    # For each connected node pair, check if there's class switching
    for i in range(nnodes):
        for j in range(nnodes):
            if sn.connmatrix[i, j] <= 0:
                continue

            # Build the class switching matrix for this node pair
            cs_matrix = np.zeros((K, K))
            has_nonzero = False

            for r in range(K):
                row_sum = 0.0
                for s in range(K):
                    src_idx = i * K + r
                    dst_idx = j * K + s
                    if src_idx < sn.rtnodes.shape[0] and dst_idx < sn.rtnodes.shape[1]:
                        prob = sn.rtnodes[src_idx, dst_idx]
                        cs_matrix[r, s] = prob
                        row_sum += prob
                        if prob > 0:
                            has_nonzero = True

                # Normalize row if it has non-zero entries
                if row_sum > 0:
                    cs_matrix[r, :] /= row_sum

            if not has_nonzero:
                continue

            # Check if matrix is not diagonal (i.e., has class switching)
            is_diagonal = True
            for r in range(K):
                for s in range(K):
                    if r != s and cs_matrix[r, s] > 1e-10:
                        is_diagonal = False
                        break
                if not is_diagonal:
                    break

            if not is_diagonal:
                cs_nodes[(i, j)] = cs_matrix

    return cs_nodes


def _write_auto_classswitch_node(
    sim: ET.Element,
    cs_name: str,
    cs_matrix: np.ndarray,
    dest_node: str,
    classnames: List[str]
) -> None:
    """
    Write an auto-generated ClassSwitch node for handling class switching in routing.

    Args:
        sim: Parent XML element (sim)
        cs_name: Name for the ClassSwitch node (e.g., "CS_Queue1_to_Delay")
        cs_matrix: K x K class switching probability matrix
        dest_node: Name of destination node to route to
        classnames: List of class names
    """
    K = len(classnames)

    node_elem = ET.SubElement(sim, 'node')
    node_elem.set('name', cs_name)

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ClassSwitch section
    cs_section = ET.SubElement(node_elem, 'section')
    cs_section.set('className', 'ClassSwitch')

    matrix_param = ET.SubElement(cs_section, 'parameter')
    matrix_param.set('array', 'true')
    matrix_param.set('classPath', 'java.lang.Object')
    matrix_param.set('name', 'matrix')

    for r in range(K):
        ref_class = ET.SubElement(matrix_param, 'refClass')
        ref_class.text = classnames[r]

        row_param = ET.SubElement(matrix_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Float')
        row_param.set('name', 'row')

        for s in range(K):
            ref_class_col = ET.SubElement(row_param, 'refClass')
            ref_class_col.text = classnames[s]

            cell_param = ET.SubElement(row_param, 'subParameter')
            cell_param.set('classPath', 'java.lang.Float')
            cell_param.set('name', 'cell')
            cell_value = ET.SubElement(cell_param, 'value')
            cell_value.text = f'{cs_matrix[r, s]:.12f}'

    # 3. Router section - always route to destination with Random strategy
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')

    routing_param = ET.SubElement(router, 'parameter')
    routing_param.set('array', 'true')
    routing_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    routing_param.set('name', 'RoutingStrategy')

    for r in range(K):
        ref_class = ET.SubElement(routing_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(routing_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
        sub_param.set('name', 'Random')


def _write_jsim_file(sn: NetworkStruct, model_path: str, options: SolverJMTOptions) -> None:
    """
    Write the network model to JSIM XML format.

    This is a simplified version supporting basic queueing networks.
    """
    M = sn.nstations
    K = sn.nclasses

    from datetime import datetime
    timestamp = datetime.now().strftime('%a %b %d %H:%M:%S %Y')

    # Create archive root element
    archive = ET.Element('archive')
    archive.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    archive.set('name', os.path.basename(model_path))
    archive.set('timestamp', timestamp)
    archive.set('xsi:noNamespaceSchemaLocation', 'Archive.xsd')

    # Create sim element inside archive
    sim = ET.SubElement(archive, 'sim')
    sim.set('disableStatisticStop', 'true')
    sim.set('logDecimalSeparator', '.')
    sim.set('logDelimiter', ';')
    sim.set('logPath', os.path.dirname(model_path))
    sim.set('logReplaceMode', '0')
    sim.set('maxEvents', '-1')
    sim.set('maxSamples', str(options.samples))
    sim.set('name', os.path.basename(model_path))
    sim.set('polling', '1.0')
    sim.set('seed', str(options.seed))
    sim.set('xsi:noNamespaceSchemaLocation', 'SIMmodeldefinition.xsd')

    # Create class definitions
    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
    classnames = sn.classnames if sn.classnames else [f'Class{i+1}' for i in range(K)]
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]
    refstat = sn.refstat.flatten() if hasattr(sn, 'refstat') and sn.refstat is not None else np.zeros(K, dtype=int)

    for r in range(K):
        class_elem = ET.SubElement(sim, 'userClass')
        class_elem.set('name', classnames[r])
        class_elem.set('priority', '0')

        if np.isinf(njobs[r]):
            # Open class - reference source is Source node
            class_elem.set('referenceSource', 'Source')
            class_elem.set('type', 'open')
        else:
            # Closed class - reference source is the station where jobs start
            ref_idx = int(refstat[r])
            ref_name = nodenames[ref_idx] if ref_idx < len(nodenames) else classnames[r] + '_RefStation'
            class_elem.set('referenceSource', ref_name)
            class_elem.set('type', 'closed')
            class_elem.set('customers', str(int(njobs[r])))

    # Detect class switches in routing - JMT requires ClassSwitch nodes for class transitions
    cs_nodes = _detect_class_switches(sn)
    # Map (source, dest) -> ClassSwitch node name
    cs_node_names = {}
    for (src_idx, dst_idx), cs_matrix in cs_nodes.items():
        cs_name = f"CS_{nodenames[src_idx]}_to_{nodenames[dst_idx]}"
        cs_node_names[(src_idx, dst_idx)] = cs_name

    # Create nodes
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    for i in range(sn.nnodes):
        node_type = sn.nodetype[i] if sn.nodetype is not None and len(sn.nodetype) > i else NodeType.QUEUE
        node_name = nodenames[i]

        node_elem = ET.SubElement(sim, 'node')
        node_elem.set('name', node_name)

        if node_type == NodeType.SOURCE:
            _write_source_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.SINK:
            _write_sink_node(node_elem, sn, classnames)
        elif node_type == NodeType.DELAY:
            _write_delay_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.QUEUE:
            _write_queue_node(node_elem, i, sn, classnames, options, cs_node_names)
        elif node_type == NodeType.ROUTER:
            _write_router_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.FORK:
            _write_fork_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.JOIN:
            _write_join_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.CLASSSWITCH:
            _write_classswitch_node(node_elem, i, sn, classnames)
        elif node_type == NodeType.PLACE:
            _write_place_node(node_elem, i, sn, classnames)
        elif node_type == NodeType.TRANSITION:
            _write_transition_node(node_elem, i, sn, classnames)

    # Create auto-generated ClassSwitch nodes for class switching in routing
    for (src_idx, dst_idx), cs_matrix in cs_nodes.items():
        cs_name = cs_node_names[(src_idx, dst_idx)]
        dest_name = nodenames[dst_idx]
        _write_auto_classswitch_node(sim, cs_name, cs_matrix, dest_name, classnames)

    # Metrics (must come before connections per JMT schema)
    for i in range(M):
        node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
        node_name = nodenames[node_idx]
        node_type = sn.nodetype[node_idx] if sn.nodetype is not None and len(sn.nodetype) > node_idx else NodeType.QUEUE

        if node_type == NodeType.SOURCE or node_type == NodeType.SINK:
            continue

        # Format alpha to avoid floating-point precision issues (e.g., 1-0.99 = 0.010000000000000009)
        alpha_str = f'{round(1 - options.conf_int, 10)}'

        for r in range(K):
            # Queue length
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'{node_name}_{classnames[r]}_QLen')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Number of Customers')
            metric.set('verbose', 'false')

            # Response time
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'{node_name}_{classnames[r]}_RespT')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Response Time')
            metric.set('verbose', 'false')

            # Utilization
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'{node_name}_{classnames[r]}_Util')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Utilization')
            metric.set('verbose', 'false')

            # Throughput
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'{node_name}_{classnames[r]}_Tput')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Throughput')
            metric.set('verbose', 'false')

    # Create connections (must come after metrics per JMT schema)
    # When class switching exists between i and j, route through ClassSwitch node:
    # i -> CS_i_to_j -> j instead of i -> j
    if sn.connmatrix is not None:
        for i in range(sn.nnodes):
            for j in range(sn.nnodes):
                if sn.connmatrix[i, j] > 0:
                    if (i, j) in cs_node_names:
                        # Route through ClassSwitch node
                        cs_name = cs_node_names[(i, j)]
                        # Connection: source -> ClassSwitch
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', nodenames[i])
                        conn.set('target', cs_name)
                        # Connection: ClassSwitch -> dest
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', cs_name)
                        conn.set('target', nodenames[j])
                    else:
                        # Direct connection (no class switching)
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', nodenames[i])
                        conn.set('target', nodenames[j])

    # Create preload section for closed classes (MATLAB writeJSIM.m lines 145-185)
    # This is essential for closed networks - without it, JMT has no jobs to simulate
    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
    has_closed_classes = any(np.isfinite(njobs[r]) and njobs[r] > 0 for r in range(K))

    if has_closed_classes:
        preload = ET.SubElement(sim, 'preload')

        # For each station (excluding Source and Join nodes)
        for ist in range(M):
            node_idx = int(sn.stationToNode[ist]) if sn.stationToNode is not None else ist
            node_type = int(sn.nodetype[node_idx]) if sn.nodetype is not None else -1

            # Skip Source (0) and Join (5) nodes
            if node_type == 0 or node_type == 5:
                continue

            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

            # Get initial state for this station
            # For closed classes, jobs start at their reference station
            station_has_jobs = False
            class_populations = []

            for r in range(K):
                if np.isfinite(njobs[r]) and njobs[r] > 0:
                    # Check if this station is the reference station for class r
                    ref_idx = int(refstat[r]) if r < len(refstat) else 0
                    if ist == ref_idx:
                        # All jobs of this class start at reference station
                        class_populations.append((classnames[r], int(njobs[r])))
                        station_has_jobs = True

            if station_has_jobs:
                station_pop = ET.SubElement(preload, 'stationPopulations')
                station_pop.set('stationName', node_name)

                for class_name, pop in class_populations:
                    class_pop = ET.SubElement(station_pop, 'classPopulation')
                    class_pop.set('population', str(pop))
                    class_pop.set('refClass', class_name)

    # Write to file
    xml_str = ET.tostring(archive, encoding='unicode')
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')

    with open(model_path, 'w') as f:
        f.write(pretty_xml)


def _write_source_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                       cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write source node section."""
    doc = node_elem

    section = ET.SubElement(node_elem, 'section')
    section.set('className', 'RandomSource')

    param = ET.SubElement(section, 'parameter')
    param.set('array', 'true')
    param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    param.set('name', 'ServiceStrategy')

    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0
    K = sn.nclasses

    for r in range(K):
        ref_class = ET.SubElement(param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
        if not np.isinf(njobs[r]):
            # Closed class - no arrivals at source
            value = ET.SubElement(sub_param, 'value')
            value.text = 'null'
        else:
            # Open class - exponential arrivals
            rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
            if np.isnan(rate) or rate <= 0:
                value = ET.SubElement(sub_param, 'value')
                value.text = 'null'
            else:
                mean = 1.0 / rate

                distr = ET.SubElement(sub_param, 'subParameter')
                distr.set('classPath', 'jmt.engine.random.Exponential')
                distr.set('name', 'Exponential')

                distr_par = ET.SubElement(sub_param, 'subParameter')
                distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
                distr_par.set('name', 'distrPar')

                lambda_param = ET.SubElement(distr_par, 'subParameter')
                lambda_param.set('classPath', 'java.lang.Double')
                lambda_param.set('name', 'lambda')
                value = ET.SubElement(lambda_param, 'value')
                value.text = str(rate)

    # ServiceTunnel and Router sections
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_sink_node(node_elem: ET.Element, sn: NetworkStruct, classnames: List[str]):
    """Write sink node section."""
    section = ET.SubElement(node_elem, 'section')
    section.set('className', 'JobSink')


def _write_delay_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                      cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write delay station node.

    Delay nodes (infinite servers) in JMT need:
    1. Queue section (input buffer)
    2. Delay section (service strategy)
    3. Router section (routing)
    """
    K = sn.nclasses
    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0

    # Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'waiting queue'  # Match Java format

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # Delay section (service)
    server = ET.SubElement(node_elem, 'section')
    server.set('className', 'Delay')

    service_param = ET.SubElement(server, 'parameter')
    service_param.set('array', 'true')
    service_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    service_param.set('name', 'ServiceStrategy')

    for r in range(K):
        ref_class = ET.SubElement(service_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(service_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
        if np.isnan(rate) or rate < 0:
            # Disabled service - null
            value = ET.SubElement(sub_param, 'value')
            value.text = 'null'
        elif rate == 0:
            # Immediate service (zero service time)
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')
        else:
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Exponential')
            distr.set('name', 'Exponential')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
            distr_par.set('name', 'distrPar')

            lambda_param = ET.SubElement(distr_par, 'subParameter')
            lambda_param.set('classPath', 'java.lang.Double')
            lambda_param.set('name', 'lambda')
            value = ET.SubElement(lambda_param, 'value')
            value.text = str(rate)

    # Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_queue_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                      options: SolverJMTOptions, cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write queue station node."""
    K = sn.nclasses
    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0

    # Get scheduling strategy
    sched = SchedStrategy.FCFS
    if sn.sched and ist in sn.sched:
        sched = sn.sched[ist]

    # Queue section
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    # 1. Size parameter
    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # 2. Drop strategies (required)
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # 3. Queue get strategy
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', _get_sched_strategy_class(sched))
    strategy_param.set('name', 'FCFSstrategy')

    # 4. Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        # TailStrategy for FCFS (default)
        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # Server section
    server = ET.SubElement(node_elem, 'section')
    # Use PSServer for processor sharing variants (PS, DPS, GPS)
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS):
        server.set('className', 'PSServer')
    else:
        server.set('className', 'Server')

    # Number of servers
    nservers = 1
    if sn.nservers is not None:
        nservers_val = sn.nservers[ist] if len(sn.nservers.shape) == 1 else sn.nservers[ist, 0]
        if np.isinf(nservers_val):
            nservers = 1000000  # Very large number for infinite servers
        else:
            nservers = int(nservers_val)

    servers_param = ET.SubElement(server, 'parameter')
    servers_param.set('classPath', 'java.lang.Integer')
    servers_param.set('name', 'maxJobs')
    value = ET.SubElement(servers_param, 'value')
    value.text = str(nservers)

    # Number of visits
    visits_param = ET.SubElement(server, 'parameter')
    visits_param.set('array', 'true')
    visits_param.set('classPath', 'java.lang.Integer')
    visits_param.set('name', 'numberOfVisits')

    for r in range(K):
        ref_class = ET.SubElement(visits_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(visits_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'numberOfVisits')
        value = ET.SubElement(sub_param, 'value')
        value.text = '1'

    # Service strategy
    service_param = ET.SubElement(server, 'parameter')
    service_param.set('array', 'true')
    service_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    service_param.set('name', 'ServiceStrategy')

    for r in range(K):
        ref_class = ET.SubElement(service_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(service_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0

        # Get process type
        procid = None
        if hasattr(sn, 'procid') and sn.procid is not None:
            if ist < sn.procid.shape[0] and r < sn.procid.shape[1]:
                procid = sn.procid[ist, r]

        if np.isnan(rate) or rate < 0:
            # Disabled service - null
            value = ET.SubElement(sub_param, 'value')
            value.text = 'null'
        elif rate == 0:
            # Immediate service (zero service time)
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')
        elif procid in (ProcessType.MAP, ProcessType.MMPP2):
            # MAP/MMPP2 distribution
            _write_map_service_distribution(sub_param, sn, ist, r)
        else:
            # Default: Exponential distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Exponential')
            distr.set('name', 'Exponential')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
            distr_par.set('name', 'distrPar')

            lambda_param = ET.SubElement(distr_par, 'subParameter')
            lambda_param.set('classPath', 'java.lang.Double')
            lambda_param.set('name', 'lambda')
            value = ET.SubElement(lambda_param, 'value')
            value.text = str(rate)

    # PSStrategy (for PS/DPS/GPS scheduling)
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS):
        ps_strategy_param = ET.SubElement(server, 'parameter')
        ps_strategy_param.set('array', 'true')
        ps_strategy_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategy')
        ps_strategy_param.set('name', 'PSStrategy')

        for r in range(K):
            ref_class = ET.SubElement(ps_strategy_param, 'refClass')
            ref_class.text = classnames[r]

            sub_param = ET.SubElement(ps_strategy_param, 'subParameter')
            if sched == SchedStrategy.PS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy')
                sub_param.set('name', 'EPSStrategy')
            elif sched == SchedStrategy.DPS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategy')
                sub_param.set('name', 'DPSStrategy')
            elif sched == SchedStrategy.GPS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategy')
                sub_param.set('name', 'GPSStrategy')

    # Service weights (required for PSServer - PS/DPS/GPS scheduling)
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS):
        weights_param = ET.SubElement(server, 'parameter')
        weights_param.set('array', 'true')
        weights_param.set('classPath', 'java.lang.Double')
        weights_param.set('name', 'serviceWeights')

        for r in range(K):
            ref_class = ET.SubElement(weights_param, 'refClass')
            ref_class.text = classnames[r]

            sub_param = ET.SubElement(weights_param, 'subParameter')
            sub_param.set('classPath', 'java.lang.Double')
            sub_param.set('name', 'serviceWeight')
            value = ET.SubElement(sub_param, 'value')
            # Get weight from schedparam (for DPS/GPS), default to 1 for PS
            weight = 1.0
            if sched in (SchedStrategy.DPS, SchedStrategy.GPS):
                if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                    if ist < sn.schedparam.shape[0] and r < sn.schedparam.shape[1]:
                        w = sn.schedparam[ist, r]
                        if not np.isnan(w) and w > 0:
                            weight = w
            value.text = str(weight)

    # Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_router_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                       cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write router node section.

    Router nodes in JMT need:
    1. Queue (input section - buffer)
    2. ServiceTunnel (middle section - pass-through)
    3. Router (output section with routing strategy)
    """
    K = sn.nclasses

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ServiceTunnel section (pass-through)
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_classswitch_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write ClassSwitch node section.

    ClassSwitch nodes in JMT need:
    1. Queue (input section - buffer)
    2. ClassSwitch (middle section - class switching matrix)
    3. Router (output section with routing strategy)
    """
    K = sn.nclasses

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ClassSwitch section (class switching matrix)
    cs_section = ET.SubElement(node_elem, 'section')
    cs_section.set('className', 'ClassSwitch')

    # Build the class switching matrix
    matrix_param = ET.SubElement(cs_section, 'parameter')
    matrix_param.set('array', 'true')
    matrix_param.set('classPath', 'java.lang.Object')
    matrix_param.set('name', 'matrix')

    # Get connections from this node to determine destination nodes
    if sn.connmatrix is not None:
        conn_i = sn.connmatrix[node_idx, :]
        jset = np.where(conn_i > 0)[0]
    else:
        jset = np.array([])

    for r in range(K):
        ref_class = ET.SubElement(matrix_param, 'refClass')
        ref_class.text = classnames[r]

        row_param = ET.SubElement(matrix_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Float')
        row_param.set('name', 'row')

        for s in range(K):
            ref_class_col = ET.SubElement(row_param, 'refClass')
            ref_class_col.text = classnames[s]

            cell_param = ET.SubElement(row_param, 'subParameter')
            cell_param.set('classPath', 'java.lang.Float')
            cell_param.set('name', 'cell')
            cell_value = ET.SubElement(cell_param, 'value')

            # Calculate class switching probability from rtnodes
            val = 0.0
            if sn.rtnodes is not None and len(jset) > 0:
                for j in jset:
                    src_idx = node_idx * K + r
                    dst_idx = int(j) * K + s
                    if src_idx < sn.rtnodes.shape[0] and dst_idx < sn.rtnodes.shape[1]:
                        val += sn.rtnodes[src_idx, dst_idx]
            elif r == s:
                # Default: keep same class
                val = 1.0
            cell_value.text = f'{val:.12f}'

    # 3. Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames)


def _write_place_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write Place node section for Petri nets.

    Place nodes in JMT use a Storage section.
    """
    K = sn.nclasses

    # Storage section
    storage = ET.SubElement(node_elem, 'section')
    storage.set('className', 'Storage')

    # Total capacity
    capacity_param = ET.SubElement(storage, 'parameter')
    capacity_param.set('classPath', 'java.lang.Integer')
    capacity_param.set('name', 'totalCapacity')
    value = ET.SubElement(capacity_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Place capacities (per-class)
    place_cap = ET.SubElement(storage, 'parameter')
    place_cap.set('array', 'true')
    place_cap.set('classPath', 'java.lang.Integer')
    place_cap.set('name', 'capacities')

    for r in range(K):
        ref_class = ET.SubElement(place_cap, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(place_cap, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'capacity')
        value = ET.SubElement(sub_param, 'value')
        value.text = '-1'  # Infinite capacity

    # Drop rules
    drop_rules = ET.SubElement(storage, 'parameter')
    drop_rules.set('array', 'true')
    drop_rules.set('classPath', 'java.lang.String')
    drop_rules.set('name', 'dropRules')

    for r in range(K):
        ref_class = ET.SubElement(drop_rules, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_rules, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropRule')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Get strategy (FCFS)
    get_strategy = ET.SubElement(storage, 'parameter')
    get_strategy.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    get_strategy.set('name', 'FCFSstrategy')

    # Put strategies
    put_strategy = ET.SubElement(storage, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'putStrategies')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ServiceTunnel section (pass-through - tokens don't need service)
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Router section (for routing to connected transitions)
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames)


def _write_transition_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write Transition node section for Petri nets.

    Transition nodes in JMT have three sections:
    1. Enabling (enabling and inhibiting conditions)
    2. Timing (mode names, servers, timing strategies)
    3. Firing (firing outcomes)
    """
    K = sn.nclasses
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    # Get transition parameters from nodeparam
    trans_param = None
    if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
        trans_param = sn.nodeparam.get(node_idx)

    # Defaults if no nodeparam
    nmodes = 1
    modenames = ['Mode1']
    firing_prio = [1]
    fire_weight = [1.0]
    nmodeservers = [1]
    enabling = [np.zeros((sn.nnodes, K))]
    inhibiting = [np.full((sn.nnodes, K), np.inf)]
    firing = [np.zeros((sn.nnodes, K))]

    if trans_param is not None:
        nmodes = getattr(trans_param, 'nmodes', 1)
        modenames = getattr(trans_param, 'modenames', ['Mode1'])
        firing_prio = getattr(trans_param, 'firingprio', [1])
        fire_weight = getattr(trans_param, 'fireweight', [1.0])
        nmodeservers = getattr(trans_param, 'nmodeservers', np.array([1]))
        enabling = getattr(trans_param, 'enabling', enabling)
        inhibiting = getattr(trans_param, 'inhibiting', inhibiting)
        firing = getattr(trans_param, 'firing', firing)

    # 1. Enabling section
    enabling_section = ET.SubElement(node_elem, 'section')
    enabling_section.set('className', 'Enabling')

    # Enabling conditions - uses TransitionMatrix structure
    enabling_param = ET.SubElement(enabling_section, 'parameter')
    enabling_param.set('array', 'true')
    enabling_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    enabling_param.set('name', 'enablingConditions')

    for m in range(nmodes):
        # Each mode has a TransitionMatrix
        mode_matrix = ET.SubElement(enabling_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'enablingCondition')

        # enablingVectors array
        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'enablingVectors')

        for k in range(sn.nnodes):
            if sn.nodetype[k] != NodeType.PLACE:
                continue

            # Check if this place has relevant entries
            has_relevant = False
            for r in range(K):
                en_val = enabling[m][k, r] if m < len(enabling) else 0
                in_val = inhibiting[m][k, r] if m < len(inhibiting) else np.inf
                if (not np.isinf(en_val) and en_val > 0) or (not np.isinf(in_val) and in_val > 0):
                    has_relevant = True
                    break

            if not has_relevant:
                continue

            # Create vector for this place
            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'enablingVector')

            # Station name
            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            # Enabling entries array
            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'enablingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'enablingEntry')
                val = ET.SubElement(entry, 'value')
                en_val = enabling[m][k, r] if m < len(enabling) else 0
                val.text = '-1' if np.isinf(en_val) else str(int(en_val))

    # Inhibiting conditions - same structure as enabling
    inhibiting_param = ET.SubElement(enabling_section, 'parameter')
    inhibiting_param.set('array', 'true')
    inhibiting_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    inhibiting_param.set('name', 'inhibitingConditions')

    for m in range(nmodes):
        mode_matrix = ET.SubElement(inhibiting_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'inhibitingCondition')

        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'inhibitingVectors')

        for k in range(sn.nnodes):
            if sn.nodetype[k] != NodeType.PLACE:
                continue

            has_relevant = False
            for r in range(K):
                in_val = inhibiting[m][k, r] if m < len(inhibiting) else np.inf
                if not np.isinf(in_val):
                    has_relevant = True
                    break

            if not has_relevant:
                continue

            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'inhibitingVector')

            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'inhibitingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'inhibitingEntry')
                val = ET.SubElement(entry, 'value')
                in_val = inhibiting[m][k, r] if m < len(inhibiting) else np.inf
                val.text = '-1' if np.isinf(in_val) else str(int(in_val))

    # 2. Timing section
    timing_section = ET.SubElement(node_elem, 'section')
    timing_section.set('className', 'Timing')

    # Mode names
    modenames_param = ET.SubElement(timing_section, 'parameter')
    modenames_param.set('array', 'true')
    modenames_param.set('classPath', 'java.lang.String')
    modenames_param.set('name', 'modeNames')

    for m in range(nmodes):
        sub_param = ET.SubElement(modenames_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'modeName')
        value = ET.SubElement(sub_param, 'value')
        value.text = modenames[m]

    # Number of servers
    servers_param = ET.SubElement(timing_section, 'parameter')
    servers_param.set('array', 'true')
    servers_param.set('classPath', 'java.lang.Integer')
    servers_param.set('name', 'numbersOfServers')

    for m in range(nmodes):
        sub_param = ET.SubElement(servers_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'numberOfServers')
        value = ET.SubElement(sub_param, 'value')
        nservers_raw = nmodeservers[m] if m < len(nmodeservers) else 1
        if np.isinf(nservers_raw) or nservers_raw >= 1000000:
            value.text = '-1'  # JMT uses -1 for infinite servers
        else:
            value.text = str(int(nservers_raw))

    # Timing strategies
    timing_strat_param = ET.SubElement(timing_section, 'parameter')
    timing_strat_param.set('array', 'true')
    timing_strat_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    timing_strat_param.set('name', 'timingStrategies')

    for m in range(nmodes):
        dist = None
        if trans_param is not None and hasattr(trans_param, 'distributions'):
            dists = trans_param.distributions
            if m < len(dists) and dists[m] is not None:
                dist = dists[m]

        if dist is not None:
            sub_param = ET.SubElement(timing_strat_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
            sub_param.set('name', 'timingStrategy')
            _write_distribution_param(sub_param, dist)
        else:
            sub_param = ET.SubElement(timing_strat_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')

    # Firing priorities
    prio_param = ET.SubElement(timing_section, 'parameter')
    prio_param.set('array', 'true')
    prio_param.set('classPath', 'java.lang.Integer')
    prio_param.set('name', 'firingPriorities')

    for m in range(nmodes):
        sub_param = ET.SubElement(prio_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'firingPriority')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(int(firing_prio[m]) if m < len(firing_prio) else 1)

    # Firing weights
    weight_param = ET.SubElement(timing_section, 'parameter')
    weight_param.set('array', 'true')
    weight_param.set('classPath', 'java.lang.Double')
    weight_param.set('name', 'firingWeights')

    for m in range(nmodes):
        sub_param = ET.SubElement(weight_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Double')
        sub_param.set('name', 'firingWeight')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(float(fire_weight[m]) if m < len(fire_weight) else 1.0)

    # 3. Firing section
    firing_section = ET.SubElement(node_elem, 'section')
    firing_section.set('className', 'Firing')

    # Firing outcomes - same structure as enabling
    firing_param = ET.SubElement(firing_section, 'parameter')
    firing_param.set('array', 'true')
    firing_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    firing_param.set('name', 'firingOutcomes')

    for m in range(nmodes):
        mode_matrix = ET.SubElement(firing_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'firingOutcome')

        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'firingVectors')

        for k in range(sn.nnodes):
            if sn.nodetype[k] != NodeType.PLACE:
                continue

            has_relevant = False
            for r in range(K):
                fire_val = firing[m][k, r] if m < len(firing) else 0
                if fire_val != 0:
                    has_relevant = True
                    break

            if not has_relevant:
                continue

            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'firingVector')

            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'firingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'firingEntry')
                val = ET.SubElement(entry, 'value')
                fire_val = firing[m][k, r] if m < len(firing) else 0
                val.text = str(int(fire_val))


def _write_distribution_param(parent: ET.Element, dist) -> None:
    """Write distribution parameters for timing strategy."""
    dist_name = dist._name if hasattr(dist, '_name') else type(dist).__name__

    if dist_name == 'Exp':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        rate = dist.get_rate() if hasattr(dist, 'get_rate') else 1.0
        value.text = str(rate)
    elif dist_name == 'Erlang':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Erlang')
        distr.set('name', 'Erlang')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
        distr_par.set('name', 'distrPar')

        alpha_param = ET.SubElement(distr_par, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        rate = dist.get_rate() if hasattr(dist, 'get_rate') else 1.0
        value.text = str(rate)

        r_param = ET.SubElement(distr_par, 'subParameter')
        r_param.set('classPath', 'java.lang.Long')
        r_param.set('name', 'r')
        value = ET.SubElement(r_param, 'value')
        try:
            order = dist.get_number_of_phases() if hasattr(dist, 'get_number_of_phases') else 1
        except NotImplementedError:
            order = 1
        value.text = str(order)
    elif dist_name == 'HyperExp':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.HyperExp')
        distr.set('name', 'Hyperexponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.HyperExpPar')
        distr_par.set('name', 'distrPar')

        # Get parameters
        p = 0.5
        lambda1 = 1.0
        lambda2 = 1.0
        if hasattr(dist, '_p') and hasattr(dist, '_lambda1') and hasattr(dist, '_lambda2'):
            p = dist._p
            lambda1 = dist._lambda1
            lambda2 = dist._lambda2

        p_param = ET.SubElement(distr_par, 'subParameter')
        p_param.set('classPath', 'java.lang.Double')
        p_param.set('name', 'p')
        value = ET.SubElement(p_param, 'value')
        value.text = str(p)

        l1_param = ET.SubElement(distr_par, 'subParameter')
        l1_param.set('classPath', 'java.lang.Double')
        l1_param.set('name', 'lambda1')
        value = ET.SubElement(l1_param, 'value')
        value.text = str(lambda1)

        l2_param = ET.SubElement(distr_par, 'subParameter')
        l2_param.set('classPath', 'java.lang.Double')
        l2_param.set('name', 'lambda2')
        value = ET.SubElement(l2_param, 'value')
        value.text = str(lambda2)
    else:
        # Default to exponential with rate 1
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        value.text = '1.0'


def _write_map_service_distribution(parent: ET.Element, sn: NetworkStruct, ist: int, r: int) -> None:
    """
    Write MAP/MMPP2 service distribution to JMT XML.

    This writes the MAPDistr format used by JMT for Markov Arrival Processes
    and Markov Modulated Poisson Processes.

    Args:
        parent: Parent XML element (serviceTimeStrategyNode)
        sn: NetworkStruct containing proc field with D0, D1 matrices
        ist: Station index
        r: Class index
    """
    # Get D0 and D1 matrices from proc
    D0 = None
    D1 = None

    if hasattr(sn, 'proc') and sn.proc is not None:
        try:
            proc = sn.proc[ist][r]
            if proc is not None and isinstance(proc, (list, tuple)) and len(proc) >= 2:
                D0 = np.asarray(proc[0], dtype=np.float64)
                D1 = np.asarray(proc[1], dtype=np.float64)
        except (IndexError, TypeError):
            pass

    if D0 is None or D1 is None:
        # Fallback to exponential if matrices not available
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
        value.text = str(rate)
        return

    # Get number of phases
    n_phases = D0.shape[0]

    # Write distribution element
    distr = ET.SubElement(parent, 'subParameter')
    distr.set('classPath', 'jmt.engine.random.MAPDistr')
    distr.set('name', 'Burst (MAP)')

    # Write parameter element
    distr_par = ET.SubElement(parent, 'subParameter')
    distr_par.set('classPath', 'jmt.engine.random.MAPPar')
    distr_par.set('name', 'distrPar')

    # Write D0 matrix
    d0_param = ET.SubElement(distr_par, 'subParameter')
    d0_param.set('array', 'true')
    d0_param.set('classPath', 'java.lang.Object')
    d0_param.set('name', 'D0')

    for k in range(n_phases):
        row_param = ET.SubElement(d0_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Object')
        row_param.set('name', 'vector')

        for j in range(n_phases):
            entry_param = ET.SubElement(row_param, 'subParameter')
            entry_param.set('classPath', 'java.lang.Double')
            entry_param.set('name', 'entry')
            value = ET.SubElement(entry_param, 'value')
            value.text = f'{D0[k, j]:.12f}'

    # Write D1 matrix
    d1_param = ET.SubElement(distr_par, 'subParameter')
    d1_param.set('array', 'true')
    d1_param.set('classPath', 'java.lang.Object')
    d1_param.set('name', 'D1')

    for k in range(n_phases):
        row_param = ET.SubElement(d1_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Object')
        row_param.set('name', 'vector')

        for j in range(n_phases):
            entry_param = ET.SubElement(row_param, 'subParameter')
            entry_param.set('classPath', 'java.lang.Double')
            entry_param.set('name', 'entry')
            value = ET.SubElement(entry_param, 'value')
            value.text = f'{D1[k, j]:.12f}'


def _write_fork_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                     cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write Fork node section.

    Fork nodes in JMT have three sections:
    1. Queue (buffer)
    2. ServiceTunnel
    3. Fork (with ForkStrategy)
    """
    K = sn.nclasses

    # Get fanOut (tasks per link), default to 1
    fan_out = 1
    if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
        if node_idx in sn.nodeparam and sn.nodeparam[node_idx] is not None:
            if isinstance(sn.nodeparam[node_idx], dict):
                fan_out = sn.nodeparam[node_idx].get('fanOut', 1)

    # 1. Queue section (buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ServiceTunnel section
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Fork section
    fork = ET.SubElement(node_elem, 'section')
    fork.set('className', 'Fork')

    # jobsPerLink parameter
    jpl_param = ET.SubElement(fork, 'parameter')
    jpl_param.set('classPath', 'java.lang.Integer')
    jpl_param.set('name', 'jobsPerLink')
    value = ET.SubElement(jpl_param, 'value')
    value.text = str(fan_out)

    # block parameter
    block_param = ET.SubElement(fork, 'parameter')
    block_param.set('classPath', 'java.lang.Integer')
    block_param.set('name', 'block')
    value = ET.SubElement(block_param, 'value')
    value.text = '-1'

    # isSimplifiedFork parameter
    simpl_param = ET.SubElement(fork, 'parameter')
    simpl_param.set('classPath', 'java.lang.Boolean')
    simpl_param.set('name', 'isSimplifiedFork')
    value = ET.SubElement(simpl_param, 'value')
    value.text = 'true'

    # ForkStrategy parameter
    strategy_param = ET.SubElement(fork, 'parameter')
    strategy_param.set('array', 'true')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.ForkStrategy')
    strategy_param.set('name', 'ForkStrategy')

    # Find outgoing connections for this fork
    outgoing_nodes = []
    if sn.connmatrix is not None:
        for j in range(sn.nnodes):
            if sn.connmatrix[node_idx, j] > 0:
                outgoing_nodes.append(j)

    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    for r in range(K):
        ref_class = ET.SubElement(strategy_param, 'refClass')
        ref_class.text = classnames[r]

        class_strat = ET.SubElement(strategy_param, 'subParameter')
        class_strat.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork')
        class_strat.set('name', 'Branch Probabilities')

        emp_array = ET.SubElement(class_strat, 'subParameter')
        emp_array.set('array', 'true')
        emp_array.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath')
        emp_array.set('name', 'EmpiricalEntryArray')

        # For each outgoing link, create an OutPath entry
        for out_node in outgoing_nodes:
            out_path = ET.SubElement(emp_array, 'subParameter')
            out_path.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath')
            out_path.set('name', 'OutPathEntry')

            # outUnitProbability
            emp_entry = ET.SubElement(out_path, 'subParameter')
            emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            emp_entry.set('name', 'outUnitProbability')

            # stationName
            station_name = ET.SubElement(emp_entry, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[out_node]

            # probability
            prob_param = ET.SubElement(emp_entry, 'subParameter')
            prob_param.set('classPath', 'java.lang.Double')
            prob_param.set('name', 'probability')
            value = ET.SubElement(prob_param, 'value')
            value.text = '1.0'

            # JobsPerLinkDis
            jpl_dis = ET.SubElement(out_path, 'subParameter')
            jpl_dis.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            jpl_dis.set('array', 'true')
            jpl_dis.set('name', 'JobsPerLinkDis')

            jpl_entry = ET.SubElement(jpl_dis, 'subParameter')
            jpl_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            jpl_entry.set('name', 'EmpiricalEntry')

            # numbers (jobs per link)
            numbers_param = ET.SubElement(jpl_entry, 'subParameter')
            numbers_param.set('classPath', 'java.lang.String')
            numbers_param.set('name', 'numbers')
            value = ET.SubElement(numbers_param, 'value')
            value.text = str(fan_out)

            # probability for this distribution
            prob_param2 = ET.SubElement(jpl_entry, 'subParameter')
            prob_param2.set('classPath', 'java.lang.Double')
            prob_param2.set('name', 'probability')
            value = ET.SubElement(prob_param2, 'value')
            value.text = '1.0'


def _write_join_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                     cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write Join node section.

    Join nodes in JMT have three sections:
    1. Join (with JoinStrategy)
    2. ServiceTunnel
    3. Router (dispatcher)
    """
    K = sn.nclasses

    # Get the number of incoming links (fanIn) - this is the number of tasks to wait for
    fan_in = 0
    if sn.connmatrix is not None:
        for i in range(sn.nnodes):
            if sn.connmatrix[i, node_idx] > 0:
                fan_in += 1

    # 1. Join section
    join = ET.SubElement(node_elem, 'section')
    join.set('className', 'Join')

    strategy_param = ET.SubElement(join, 'parameter')
    strategy_param.set('array', 'true')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.JoinStrategy')
    strategy_param.set('name', 'JoinStrategy')

    for r in range(K):
        ref_class = ET.SubElement(strategy_param, 'refClass')
        ref_class.text = classnames[r]

        # Default: Standard Join (wait for all tasks)
        join_strat = ET.SubElement(strategy_param, 'subParameter')
        join_strat.set('classPath', 'jmt.engine.NetStrategies.JoinStrategies.NormalJoin')
        join_strat.set('name', 'Standard Join')

        req_param = ET.SubElement(join_strat, 'subParameter')
        req_param.set('classPath', 'java.lang.Integer')
        req_param.set('name', 'numRequired')
        value = ET.SubElement(req_param, 'value')
        value.text = str(max(1, fan_in))  # Number of tasks to wait for

    # 2. ServiceTunnel section
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_routing_strategy(router: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                            cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write routing strategy for a node.

    Handles different routing strategies (RAND, RROBIN, JSQ, PROB, etc.) by using
    the appropriate JMT strategy class.

    When class switching exists between this node and a destination, routes to
    the ClassSwitch node instead of the direct destination.

    Args:
        router: XML element for the Router section
        node_idx: Index of the current node
        sn: NetworkStruct object
        classnames: List of class names
        cs_node_names: Optional dict mapping (src_idx, dst_idx) to ClassSwitch node names
    """
    from ...sn.network_struct import RoutingStrategy

    if cs_node_names is None:
        cs_node_names = {}

    K = sn.nclasses
    M = sn.nnodes
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(M)]

    # Check if we have routing probability data
    has_rtnodes = hasattr(sn, 'rtnodes') and sn.rtnodes is not None and sn.rtnodes.size > 0
    has_connmatrix = hasattr(sn, 'connmatrix') and sn.connmatrix is not None
    has_routing = hasattr(sn, 'routing') and sn.routing is not None

    param = ET.SubElement(router, 'parameter')
    param.set('array', 'true')
    param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    param.set('name', 'RoutingStrategy')

    for r in range(K):
        ref_class = ET.SubElement(param, 'refClass')
        ref_class.text = classnames[r]

        # Get routing strategy for this node/class
        strategy = RoutingStrategy.RAND  # Default
        if has_routing and node_idx < sn.routing.shape[0] and r < sn.routing.shape[1]:
            strategy = RoutingStrategy(int(sn.routing[node_idx, r]))

        # Handle different routing strategies
        if strategy == RoutingStrategy.RROBIN:
            # Round Robin strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy')
            sub_param.set('name', 'Round Robin')

        elif strategy == RoutingStrategy.WRROBIN:
            # Weighted Round Robin strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightedRoundRobinStrategy')
            sub_param.set('name', 'Weighted Round Robin')

            weight_array = ET.SubElement(sub_param, 'subParameter')
            weight_array.set('array', 'true')
            weight_array.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry')
            weight_array.set('name', 'WeightEntryArray')

            # Add weight entries for connected nodes
            if has_connmatrix:
                for j in range(M):
                    if sn.connmatrix[node_idx, j] > 0:
                        weight_entry = ET.SubElement(weight_array, 'subParameter')
                        weight_entry.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry')
                        weight_entry.set('name', 'WeightEntry')

                        station_param = ET.SubElement(weight_entry, 'subParameter')
                        station_param.set('classPath', 'java.lang.String')
                        station_param.set('name', 'stationName')
                        station_value = ET.SubElement(station_param, 'value')
                        station_value.text = nodenames[j]

                        weight_param = ET.SubElement(weight_entry, 'subParameter')
                        weight_param.set('classPath', 'java.lang.Integer')
                        weight_param.set('name', 'weight')
                        weight_value = ET.SubElement(weight_param, 'value')
                        # Default weight of 1 for each destination
                        weight_value.text = '1'

        elif strategy == RoutingStrategy.JSQ:
            # Join Shortest Queue strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.ShortestQueueLengthRoutingStrategy')
            sub_param.set('name', 'Join the Shortest Queue (JSQ)')

        elif strategy == RoutingStrategy.KCHOICES:
            # Power of K choices strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy')
            sub_param.set('name', 'Power of k')

            k_param = ET.SubElement(sub_param, 'subParameter')
            k_param.set('classPath', 'java.lang.Integer')
            k_param.set('name', 'k')
            k_value = ET.SubElement(k_param, 'value')
            k_value.text = '2'  # Default k=2

            mem_param = ET.SubElement(sub_param, 'subParameter')
            mem_param.set('classPath', 'java.lang.Boolean')
            mem_param.set('name', 'withMemory')
            mem_value = ET.SubElement(mem_param, 'value')
            mem_value.text = 'false'

        elif strategy == RoutingStrategy.PROB:
            # Probabilistic routing using EmpiricalStrategy
            # When class switching exists, sum probabilities across all destination classes
            # and route to ClassSwitch node instead of direct destination
            routing_probs = []
            if has_rtnodes and has_connmatrix:
                row_idx = node_idx * K + r
                if row_idx < sn.rtnodes.shape[0]:
                    for j in range(M):
                        if sn.connmatrix[node_idx, j] > 0:
                            # Check if there's a ClassSwitch node for this edge
                            if (node_idx, j) in cs_node_names:
                                # Sum probabilities across ALL destination classes
                                total_prob = 0.0
                                for s in range(K):
                                    col_idx = j * K + s
                                    if col_idx < sn.rtnodes.shape[1]:
                                        total_prob += sn.rtnodes[row_idx, col_idx]
                                if total_prob > 0:
                                    routing_probs.append((cs_node_names[(node_idx, j)], total_prob))
                            else:
                                # No class switching - use same-class probability
                                col_idx = j * K + r
                                if col_idx < sn.rtnodes.shape[1]:
                                    prob = sn.rtnodes[row_idx, col_idx]
                                    if prob > 0:
                                        routing_probs.append((nodenames[j], prob))

            if routing_probs:
                # Normalize probabilities to sum to 1.0 if needed
                total = sum(p for _, p in routing_probs)
                if total > 0 and abs(total - 1.0) > 1e-10:
                    routing_probs = [(name, prob / total) for name, prob in routing_probs]

                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy')
                sub_param.set('name', 'Probabilities')

                emp_array = ET.SubElement(sub_param, 'subParameter')
                emp_array.set('array', 'true')
                emp_array.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                emp_array.set('name', 'EmpiricalEntryArray')

                for dest_name, prob in routing_probs:
                    emp_entry = ET.SubElement(emp_array, 'subParameter')
                    emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                    emp_entry.set('name', 'EmpiricalEntry')

                    station_param = ET.SubElement(emp_entry, 'subParameter')
                    station_param.set('classPath', 'java.lang.String')
                    station_param.set('name', 'stationName')
                    station_value = ET.SubElement(station_param, 'value')
                    station_value.text = dest_name

                    prob_param = ET.SubElement(emp_entry, 'subParameter')
                    prob_param.set('classPath', 'java.lang.Double')
                    prob_param.set('name', 'probability')
                    prob_value = ET.SubElement(prob_param, 'value')
                    prob_value.text = f'{prob:.12f}'
            else:
                # Fallback to Random if no probabilities defined
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
                sub_param.set('name', 'Random')

        elif strategy == RoutingStrategy.DISABLED:
            # Disabled routing
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy')
            sub_param.set('name', 'Disabled')

        else:
            # Default: Try to use EmpiricalStrategy with explicit probabilities if available
            # This handles RAND strategy when routing probabilities are defined
            # When class switching exists, sum probabilities and route to ClassSwitch node
            routing_probs = []
            if has_rtnodes and has_connmatrix:
                row_idx = node_idx * K + r
                if row_idx < sn.rtnodes.shape[0]:
                    # Iterate over all nodes (not just stations) to include Sink
                    n_nodes = sn.nnodes if hasattr(sn, 'nnodes') else sn.connmatrix.shape[1]
                    for j in range(n_nodes):
                        if j < sn.connmatrix.shape[1] and sn.connmatrix[node_idx, j] > 0:
                            # Check if there's a ClassSwitch node for this edge
                            if (node_idx, j) in cs_node_names:
                                # Sum probabilities across ALL destination classes
                                total_prob = 0.0
                                for s in range(K):
                                    col_idx = j * K + s
                                    if col_idx < sn.rtnodes.shape[1]:
                                        total_prob += sn.rtnodes[row_idx, col_idx]
                                if total_prob > 0:
                                    routing_probs.append((cs_node_names[(node_idx, j)], total_prob))
                            else:
                                # No class switching - use same-class probability
                                col_idx = j * K + r
                                if col_idx < sn.rtnodes.shape[1]:
                                    prob = sn.rtnodes[row_idx, col_idx]
                                    if prob > 0:
                                        routing_probs.append((nodenames[j], prob))

            if routing_probs:
                # Normalize probabilities to sum to 1.0
                # This is necessary for RAND strategy where rtnodes contains connection weights (1.0)
                # rather than actual routing probabilities
                total = sum(p for _, p in routing_probs)
                if total > 0 and abs(total - 1.0) > 1e-10:
                    routing_probs = [(name, prob / total) for name, prob in routing_probs]

                # Use EmpiricalStrategy with explicit probabilities
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy')
                sub_param.set('name', 'Probabilities')

                emp_array = ET.SubElement(sub_param, 'subParameter')
                emp_array.set('array', 'true')
                emp_array.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                emp_array.set('name', 'EmpiricalEntryArray')

                for dest_name, prob in routing_probs:
                    emp_entry = ET.SubElement(emp_array, 'subParameter')
                    emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                    emp_entry.set('name', 'EmpiricalEntry')

                    station_param = ET.SubElement(emp_entry, 'subParameter')
                    station_param.set('classPath', 'java.lang.String')
                    station_param.set('name', 'stationName')
                    station_value = ET.SubElement(station_param, 'value')
                    station_value.text = dest_name

                    prob_param = ET.SubElement(emp_entry, 'subParameter')
                    prob_param.set('classPath', 'java.lang.Double')
                    prob_param.set('name', 'probability')
                    prob_value = ET.SubElement(prob_param, 'value')
                    prob_value.text = f'{prob:.12f}'
            else:
                # Fallback to Random strategy
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
                sub_param.set('name', 'Random')


def _parse_jsim_results(result_path: str, sn: NetworkStruct) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse JMT simulation results from XML output.

    Returns:
        Tuple of (Q, U, R, T) matrices
    """
    M = sn.nstations
    K = sn.nclasses

    Q = np.full((M, K), np.nan)
    U = np.full((M, K), np.nan)
    R = np.full((M, K), np.nan)
    T = np.full((M, K), np.nan)

    if not os.path.exists(result_path):
        return Q, U, R, T

    try:
        tree = ET.parse(result_path)
        root = tree.getroot()

        classnames = sn.classnames if sn.classnames else [f'Class{i+1}' for i in range(K)]
        nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

        for measure in root.iter('measure'):
            measure_type = measure.get('measureType', measure.get('type', ''))
            node_name = measure.get('station', measure.get('referenceNode', ''))
            class_name = measure.get('class', measure.get('referenceUserClass', ''))
            mean_value = measure.get('meanValue', '0')
            # Note: We accept results even when successful="false" as JMT still provides
            # valid mean values - it just means the precision target wasn't met

            # Find station index
            station_idx = -1
            for i in range(M):
                node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
                if node_idx < len(nodenames) and nodenames[node_idx] == node_name:
                    station_idx = i
                    break

            if station_idx < 0:
                continue

            # Find class index
            class_idx = -1
            for r in range(K):
                if r < len(classnames) and classnames[r] == class_name:
                    class_idx = r
                    break

            if class_idx < 0:
                continue

            try:
                value = float(mean_value)
            except ValueError:
                continue

            if 'Number of Customers' in measure_type or 'QLen' in measure_type:
                Q[station_idx, class_idx] = value
            elif 'Utilization' in measure_type or 'Util' in measure_type:
                U[station_idx, class_idx] = value
            elif 'Response Time' in measure_type or 'RespT' in measure_type:
                R[station_idx, class_idx] = value
            elif 'Throughput' in measure_type or 'Tput' in measure_type:
                T[station_idx, class_idx] = value

    except Exception as e:
        pass  # Return NaN-filled matrices on parse error

    return Q, U, R, T


def solver_jmt(
    sn: NetworkStruct,
    options: Optional[SolverJMTOptions] = None
) -> SolverJMTReturn:
    """
    JMT solver handler - calls JMT via subprocess.

    Performs discrete-event simulation using JMT by:
    1. Writing the model to JSIM XML format
    2. Calling JMT via subprocess
    3. Parsing the results

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverJMTReturn with all performance metrics

    Raises:
        RuntimeError: If JMT is not available or fails
    """
    start_time = time.time()

    if options is None:
        options = SolverJMTOptions()

    if not is_jmt_available():
        raise RuntimeError(
            "SolverJMT requires Java and JMT.jar.\n"
            "Ensure Java is installed and JMT.jar is in the common/ directory."
        )

    jmt_path = _get_jmt_jar_path()

    M = sn.nstations
    K = sn.nclasses

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix='jmt_')

    try:
        model_path = os.path.join(temp_dir, 'model.jsimg')
        result_path = model_path + '-result.jsim'  # JMT creates .jsimg-result.jsim

        # Write model to JSIM format
        _write_jsim_file(sn, model_path, options)

        # Build command
        cmd = [
            'java',
            '-cp', jmt_path,
            'jmt.commandline.Jmt',
            'sim',
            model_path,
        ]

        if options.verbose:
            print(f"SolverJMT command: {' '.join(cmd)}")

        # Execute command
        result = subprocess.run(
            cmd,
            capture_output=True,
            cwd=temp_dir,
            timeout=600  # 10 minute timeout
        )

        if result.returncode != 0:
            stderr = result.stderr.decode('utf-8', errors='ignore')
            raise RuntimeError(f"JMT simulation failed: {stderr}")

        # Parse results
        Q, U, R, T = _parse_jsim_results(result_path, sn)

        # Calculate arrival rates and other metrics
        A = np.zeros((M, K))
        W = R.copy()

        # System throughput (sum at reference stations)
        X = np.zeros((1, K))
        refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(K, dtype=int)
        for r in range(K):
            if int(refstat[r]) < M:
                X[0, r] = T[int(refstat[r]), r]

        # Cycle times
        C = np.zeros((1, K))
        njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
        for r in range(K):
            if not np.isinf(njobs[r]) and X[0, r] > 0:
                C[0, r] = njobs[r] / X[0, r]

        runtime = time.time() - start_time

        return SolverJMTReturn(
            Q=Q,
            U=U,
            R=R,
            T=T,
            A=A,
            W=W,
            C=C,
            X=X,
            runtime=runtime,
            method='jsim'
        )

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
