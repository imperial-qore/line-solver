"""
Model format converters for LINE.

This module provides functions to convert between different model representations,
including NetworkStruct to Network conversions.

Port from:
    - matlab/src/io/QN2LINE.m
    - matlab/src/io/LQN2QN.m
    - matlab/src/io/QN2LQN.m
"""

import numpy as np
from typing import Any, Optional, Dict, List, Tuple
from dataclasses import dataclass

from ..mam import map_mean, map_scv


def qn2line(sn: Any, model_name: str = 'model') -> Dict[str, Any]:
    """
    Convert NetworkStruct (QN) to a Network model representation.

    Creates a Network model from a NetworkStruct, reconstructing all nodes,
    classes, service processes, and routing.

    Args:
        sn: NetworkStruct object (from getStruct())
        model_name: Name for the created model

    Returns:
        Dictionary containing network model specification that can be used
        to construct a Network object.

    References:
        MATLAB: matlab/src/io/QN2LINE.m
    """
    M = sn.nstations  # number of stations
    K = sn.nclasses   # number of classes
    rt = sn.rt if hasattr(sn, 'rt') else None
    NK = sn.njobs if hasattr(sn, 'njobs') else np.zeros(K)
    Ktrue = np.count_nonzero(NK)  # classes that are not artificial

    # Result structure
    result = {
        'name': model_name,
        'nodes': [],
        'classes': [],
        'processes': [],
        'routing': {},
    }

    # Track source/sink
    has_sink = False
    id_source = None

    # Create nodes
    PH = sn.proc if hasattr(sn, 'proc') else None

    for ist in range(M):
        sched = sn.sched[ist] if hasattr(sn, 'sched') else None
        sched_name = sched.name if hasattr(sched, 'name') else str(sched)
        node_name = sn.nodenames[ist] if hasattr(sn, 'nodenames') else f'Station{ist}'

        node_spec = {
            'id': ist,
            'name': node_name,
            'scheduling': sched_name,
        }

        if sched_name == 'INF':
            node_spec['type'] = 'Delay'
        elif sched_name == 'FORK':
            node_spec['type'] = 'Fork'
        elif sched_name == 'EXT':
            node_spec['type'] = 'Source'
            id_source = ist
            has_sink = True
            # Add sink
            result['nodes'].append({
                'id': M,
                'name': 'Sink',
                'type': 'Sink',
            })
        else:
            node_spec['type'] = 'Queue'
            node_spec['servers'] = int(sn.nservers[ist]) if hasattr(sn, 'nservers') else 1

        result['nodes'].append(node_spec)

    # Create classes
    for k in range(K):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        class_spec = {
            'id': k,
            'name': class_name,
        }

        if k < Ktrue:
            if np.isinf(NK[k]):
                class_spec['type'] = 'open'
                class_spec['population'] = np.inf
            else:
                class_spec['type'] = 'closed'
                class_spec['population'] = int(NK[k])
                class_spec['refstation'] = int(sn.refstat[k]) if hasattr(sn, 'refstat') else 0
        else:
            # Artificial class - find first station with non-null rate
            iref = 0
            if PH is not None:
                for ist in range(M):
                    if ist < len(PH) and PH[ist] is not None:
                        if k < len(PH[ist]) and PH[ist][k] is not None:
                            if hasattr(PH[ist][k], '__getitem__') and len(PH[ist][k]) > 0:
                                if np.sum(np.abs(PH[ist][k][0])) > 0:
                                    iref = ist
                                    break

            if np.isinf(NK[k]):
                class_spec['type'] = 'open'
                class_spec['population'] = np.inf
            else:
                class_spec['type'] = 'closed'
                class_spec['population'] = int(NK[k])
                class_spec['refstation'] = iref

        result['classes'].append(class_spec)

        # Create service/arrival processes for this class
        for ist in range(M):
            if PH is not None and ist < len(PH) and PH[ist] is not None:
                if k < len(PH[ist]) and PH[ist][k] is not None:
                    try:
                        scv_ik = map_scv(PH[ist][k])
                        mean_ik = map_mean(PH[ist][k])
                    except Exception:
                        continue

                    sched = sn.sched[ist] if hasattr(sn, 'sched') else None
                    sched_name = sched.name if hasattr(sched, 'name') else str(sched)

                    rate = sn.rates[ist, k] if hasattr(sn, 'rates') else 1.0 / mean_ik

                    process_spec = {
                        'station': ist,
                        'class': k,
                        'mean': mean_ik,
                        'scv': scv_ik,
                        'rate': rate,
                    }

                    if sched_name == 'EXT':
                        process_spec['process_type'] = 'arrival'
                        if np.isnan(rate):
                            process_spec['distribution'] = 'Disabled'
                        elif rate == 0:
                            process_spec['distribution'] = 'Immediate'
                        else:
                            process_spec['distribution'] = 'APH'
                    elif sched_name != 'FORK':
                        process_spec['process_type'] = 'service'
                        if np.isnan(rate):
                            process_spec['distribution'] = 'Disabled'
                        elif rate == 0:
                            process_spec['distribution'] = 'Immediate'
                        else:
                            process_spec['distribution'] = 'APH'

                    result['processes'].append(process_spec)

    # Create routing matrix
    if rt is not None:
        for k in range(K):
            for c in range(K):
                result['routing'][(k, c)] = np.zeros((M + (1 if has_sink else 0),
                                                       M + (1 if has_sink else 0)))
                for ist in range(M):
                    for m in range(M):
                        if has_sink and m == id_source:
                            # Direct to sink instead of source
                            result['routing'][(k, c)][ist, M] = rt[ist * K + k, m * K + c]
                        else:
                            result['routing'][(k, c)][ist, m] = rt[ist * K + k, m * K + c]

    return result


def line2qn(model: Any) -> Any:
    """
    Convert a Network model to NetworkStruct.

    This is essentially calling model.getStruct().

    Args:
        model: Network model object

    Returns:
        NetworkStruct object
    """
    if hasattr(model, 'getStruct'):
        return model.getStruct()
    return model


@dataclass
class LQNTask:
    """Layered Queueing Network task specification."""
    name: str
    entries: List[Dict[str, Any]]
    host: Optional[str] = None
    multiplicity: int = 1
    scheduling: str = 'ref'


@dataclass
class LQNProcessor:
    """Layered Queueing Network processor specification."""
    name: str
    tasks: List[LQNTask]
    multiplicity: int = 1
    scheduling: str = 'fcfs'


@dataclass
class LQNModel:
    """Layered Queueing Network model."""
    name: str
    processors: List[LQNProcessor]
    calls: List[Dict[str, Any]]


def qn2lqn(model: Any) -> LQNModel:
    """
    Convert a Queueing Network to Layered Queueing Network representation.

    Creates an LQN model from a QN model by mapping:
    - Delay stations to reference tasks
    - Queue stations to tasks with entries
    - Class routing to synchronous calls

    Args:
        model: Network model or NetworkStruct

    Returns:
        LQNModel representation

    References:
        MATLAB: matlab/src/io/QN2LQN.m
    """
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    processors = []
    calls = []

    # Create one processor per station
    for ist in range(sn.nstations):
        sched = sn.sched[ist] if hasattr(sn, 'sched') else None
        sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'
        node_name = sn.nodenames[ist] if hasattr(sn, 'nodenames') else f'Station{ist}'

        # Skip source
        if sched_name == 'EXT':
            continue

        # Create entries for each class served at this station
        entries = []
        for k in range(sn.nclasses):
            if hasattr(sn, 'rates') and not np.isnan(sn.rates[ist, k]) and sn.rates[ist, k] > 0:
                class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
                entry = {
                    'name': f'{node_name}_{class_name}',
                    'service_time': 1.0 / sn.rates[ist, k],
                    'class': k,
                }
                entries.append(entry)

        if entries:
            task = LQNTask(
                name=f'Task_{node_name}',
                entries=entries,
                host=f'Proc_{node_name}',
                multiplicity=int(sn.nservers[ist]) if hasattr(sn, 'nservers') else 1,
                scheduling='ref' if sched_name == 'INF' else 'fcfs',
            )

            processor = LQNProcessor(
                name=f'Proc_{node_name}',
                tasks=[task],
                multiplicity=1,
                scheduling='fcfs',
            )
            processors.append(processor)

    # Create calls based on routing
    if hasattr(sn, 'rt') and sn.rt is not None:
        for ist in range(sn.nstations):
            for k in range(sn.nclasses):
                for m in range(sn.nstations):
                    for c in range(sn.nclasses):
                        prob = sn.rt[ist * sn.nclasses + k, m * sn.nclasses + c]
                        if prob > 0:
                            calls.append({
                                'from_station': ist,
                                'from_class': k,
                                'to_station': m,
                                'to_class': c,
                                'probability': prob,
                                'type': 'sync',
                            })

    model_name = sn.name if hasattr(sn, 'name') else 'lqn_model'
    return LQNModel(name=model_name, processors=processors, calls=calls)


def lqn2qn(lqn_model: LQNModel) -> Dict[str, Any]:
    """
    Convert a Layered Queueing Network to Queueing Network representation.

    Creates a QN specification from an LQN model.

    Args:
        lqn_model: LQNModel object

    Returns:
        Dictionary with QN model specification

    References:
        MATLAB: matlab/src/io/LQN2QN.m
    """
    result = {
        'name': lqn_model.name,
        'nodes': [],
        'classes': [],
        'processes': [],
        'routing': {},
    }

    # Map each processor to a station
    station_id = 0
    class_id = 0
    entry_to_class = {}

    for proc in lqn_model.processors:
        for task in proc.tasks:
            node_spec = {
                'id': station_id,
                'name': task.name,
                'type': 'Delay' if task.scheduling == 'ref' else 'Queue',
                'scheduling': 'INF' if task.scheduling == 'ref' else 'FCFS',
                'servers': task.multiplicity,
            }
            result['nodes'].append(node_spec)

            for entry in task.entries:
                class_spec = {
                    'id': class_id,
                    'name': entry['name'],
                    'type': 'closed',
                    'population': 1,
                    'refstation': station_id if task.scheduling == 'ref' else 0,
                }
                result['classes'].append(class_spec)

                process_spec = {
                    'station': station_id,
                    'class': class_id,
                    'mean': entry.get('service_time', 1.0),
                    'scv': 1.0,
                    'rate': 1.0 / entry.get('service_time', 1.0),
                    'process_type': 'service',
                    'distribution': 'Exp',
                }
                result['processes'].append(process_spec)

                entry_to_class[entry['name']] = class_id
                class_id += 1

            station_id += 1

    return result


@dataclass
class RandomEnvironmentModel:
    """Random Environment model specification.

    Represents a queueing network with environment-modulated service rates.
    """
    name: str
    stages: List[Dict[str, Any]]
    transitions: List[Dict[str, Any]]


@dataclass
class MMPP2Params:
    """MMPP2 distribution parameters.

    D0: Phase transition matrix (off-diagonal elements define transitions)
    D1: Service rate matrix (diagonal elements define service rates per phase)
    """
    D0: np.ndarray
    D1: np.ndarray


def mapqn2renv(model: Any, options: Optional[Dict] = None):
    """
    Transform a queueing network with MMPP service into a random environment model.

    This function transforms a queueing network where servers use MMPP2
    (2-phase Markov Modulated Poisson Process) service distributions into
    a random environment model with exponential services. The environment
    has two stages (one per MMPP phase) with transitions defined by the
    MMPP D0 matrix.

    Args:
        model: Network with MMPP2 service distributions
        options: Optional configuration dictionary (reserved for future use)

    Returns:
        RandomEnvironmentModel with exponential services modulated by MMPP phases

    Raises:
        ValueError: If no MMPP2 distributions found or invalid parameters

    The transformation works as follows:
    - Input MMPP2 has D0 (phase transitions) and D1 (diagonal service rates)
    - Output has 2 environment stages with exponential services
    - Stage 0 uses service rate λ₀ = D1(1,1)
    - Stage 1 uses service rate λ₁ = D1(2,2)
    - Environment transitions: σ₀₁ = D0(1,2), σ₁₀ = D0(2,1)

    Example:
        >>> # Create network with MMPP2 service
        >>> model = Network('MMPP_Queue')
        >>> # ... setup with MMPP2 service distributions
        >>> env_model = mapqn2renv(model)

    References:
        MATLAB: matlab/src/io/MAPQN2RENV.m
    """
    if options is None:
        options = {}

    # Get network structure
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    # Phase 1: Validate and extract MMPP parameters
    mmpp_params = _validate_and_extract_mmpp(model, sn)

    if mmpp_params is None:
        raise ValueError('Network must contain at least one MMPP2 service distribution')

    # Phase 2: Extract MMPP Parameters
    D0 = mmpp_params.D0
    D1 = mmpp_params.D1

    # Extract service rates from D1 diagonal
    lambda0 = D1[0, 0]
    lambda1 = D1[1, 1]

    # Extract transition rates from D0 off-diagonal
    sigma01 = D0[0, 1]
    sigma10 = D0[1, 0]

    # Validate rates are non-negative
    if lambda0 < 0 or lambda1 < 0 or sigma01 < 0 or sigma10 < 0:
        raise ValueError('All extracted rates must be non-negative')

    # Phase 3: Create Environment Model
    model_name = getattr(model, 'name', 'model') if hasattr(model, 'name') else 'model'
    env_name = f'{model_name}_RENV'

    # Phase 4: Build Stage Networks
    stages = []

    # Stage 0 (Phase 0) - uses lambda0 service rate
    stage0 = _build_stage_network_spec(sn, 'Phase0', lambda0)
    stages.append(stage0)

    # Stage 1 (Phase 1) - uses lambda1 service rate
    stage1 = _build_stage_network_spec(sn, 'Phase1', lambda1)
    stages.append(stage1)

    # Phase 5: Add Environment Transitions
    transitions = []

    if sigma01 > 0:
        transitions.append({
            'from_stage': 'Phase0',
            'to_stage': 'Phase1',
            'rate': sigma01,
            'distribution': 'Exp',
        })

    if sigma10 > 0:
        transitions.append({
            'from_stage': 'Phase1',
            'to_stage': 'Phase0',
            'rate': sigma10,
            'distribution': 'Exp',
        })

    # Import Environment here to avoid circular imports
    from line_solver.environment import Environment
    from line_solver.distributions import Exp

    # Create actual Environment object
    env = Environment('MAPQN_Env', 2)

    # Build stage networks and add to environment
    stage_net0 = _build_stage_network(model, 'Phase0', lambda0)
    env.add_stage(0, 'Phase0', 'item', stage_net0)

    stage_net1 = _build_stage_network(model, 'Phase1', lambda1)
    env.add_stage(1, 'Phase1', 'item', stage_net1)

    # Add transitions
    if sigma01 > 0:
        env.add_transition(0, 1, Exp(sigma01))

    if sigma10 > 0:
        env.add_transition(1, 0, Exp(sigma10))

    return env


def _validate_and_extract_mmpp(model: Any, sn: Any) -> Optional[MMPP2Params]:
    """
    Validate MMPP2 distributions and extract parameters.

    Args:
        model: Network model
        sn: NetworkStruct

    Returns:
        MMPP2Params if MMPP2 found, None otherwise
    """
    first_mmpp2 = None

    # Check if model has nodes attribute (full Network object)
    if hasattr(model, 'get_nodes'):
        nodes = model.get_nodes()
    elif hasattr(model, 'nodes'):
        nodes = model.nodes
    elif hasattr(model, 'stations'):
        nodes = model.stations
    else:
        # Working with NetworkStruct only - check proc array for MAP distributions
        if hasattr(sn, 'proc') and sn.proc is not None:
            for ist in range(sn.nstations):
                if ist < len(sn.proc) and sn.proc[ist] is not None:
                    for k in range(sn.nclasses):
                        if k < len(sn.proc[ist]) and sn.proc[ist][k] is not None:
                            proc = sn.proc[ist][k]
                            # Check if it's a MAP (has D0 and D1 components)
                            if hasattr(proc, '__len__') and len(proc) >= 2:
                                try:
                                    D0 = np.array(proc[0])
                                    D1 = np.array(proc[1])

                                    # Check if D1 is diagonal (MMPP property)
                                    if _is_diagonal(D1):
                                        if first_mmpp2 is None:
                                            first_mmpp2 = MMPP2Params(D0=D0, D1=D1)
                                except (TypeError, ValueError):
                                    continue
        return first_mmpp2

    # Iterate through nodes to find MMPP2 distributions
    for node in nodes:
        node_class_name = node.__class__.__name__ if hasattr(node, '__class__') else ''

        # Check if node is a Queue or Delay
        if node_class_name not in ('Queue', 'Delay'):
            continue

        # Try to get service distributions - handle both native Python and wrapper styles
        distributions = []

        # Native Python style: get_service method
        classes = model.get_classes() if hasattr(model, 'get_classes') else (model.classes if hasattr(model, 'classes') else [])
        if hasattr(node, 'get_service') and classes:
            for job_class in classes:
                try:
                    dist = node.get_service(job_class)
                    if dist is not None:
                        distributions.append(dist)
                except Exception:
                    pass

        # Wrapper style: server.serviceProcess
        elif hasattr(node, 'server') and hasattr(node.server, 'serviceProcess'):
            service_processes = node.server.serviceProcess
            if service_processes:
                for service_list in service_processes:
                    if service_list is None:
                        continue
                    if hasattr(service_list, '__len__') and len(service_list) > 0:
                        dist = service_list[-1] if hasattr(service_list, '__getitem__') else service_list
                    else:
                        dist = service_list
                    distributions.append(dist)

        # Check each distribution for MMPP2
        for dist in distributions:
            if dist is None:
                continue

            dist_class = dist.__class__.__name__ if hasattr(dist, '__class__') else ''

            if 'MMPP2' in dist_class or 'MMPP' in dist_class:
                # Get D0 and D1 matrices - try property access first (native), then method call
                D0, D1 = None, None
                if hasattr(dist, 'D0') and hasattr(dist, 'D1'):
                    D0 = dist.D0
                    D1 = dist.D1
                elif hasattr(dist, 'D'):
                    D0 = dist.D(0) if callable(dist.D) else dist.D[0]
                    D1 = dist.D(1) if callable(dist.D) else dist.D[1]

                if D0 is not None and D1 is not None and first_mmpp2 is None:
                    first_mmpp2 = MMPP2Params(
                        D0=np.array(D0),
                        D1=np.array(D1),
                    )
                    return first_mmpp2
            elif 'MAP' in dist_class:
                # Check if it's an MMPP (diagonal D1)
                if hasattr(dist, 'D'):
                    D0 = dist.D(0) if callable(dist.D) else dist.D[0]
                    D1 = dist.D(1) if callable(dist.D) else dist.D[1]

                    D1_arr = np.array(D1)
                    if _is_diagonal(D1_arr):
                        if first_mmpp2 is None:
                            first_mmpp2 = MMPP2Params(
                                D0=np.array(D0),
                                D1=D1_arr,
                            )
                    else:
                        from .logging import line_warning
                        line_warning('mapqn2renv',
                                   'Generic MAP detected. Only MMPP (diagonal D1) is supported.')

    return first_mmpp2


def _is_diagonal(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is diagonal within tolerance."""
    if matrix.ndim != 2:
        return False
    n, m = matrix.shape
    if n != m:
        return False
    for i in range(n):
        for j in range(m):
            if i != j and abs(matrix[i, j]) > tol:
                return False
    return True


def _build_stage_network_spec(sn: Any, stage_name: str, exp_rate: float) -> Dict[str, Any]:
    """
    Build a stage network specification with exponential services.

    Args:
        sn: NetworkStruct from original model
        stage_name: Name for this stage
        exp_rate: Exponential service rate to use (replacing MMPP)

    Returns:
        Dictionary with stage network specification
    """
    stage_spec = {
        'name': stage_name,
        'exp_rate': exp_rate,
        'nodes': [],
        'classes': [],
        'processes': [],
    }

    # Clone node specifications
    for i in range(sn.nnodes):
        node_name = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
        node_type = sn.nodetype[i] if hasattr(sn, 'nodetype') else None
        type_name = node_type.name if hasattr(node_type, 'name') else 'QUEUE'

        node_spec = {
            'id': i,
            'name': node_name,
            'type': type_name,
        }

        if type_name == 'QUEUE':
            ist = sn.nodeToStation[i] if hasattr(sn, 'nodeToStation') else i
            node_spec['servers'] = int(sn.nservers[ist]) if hasattr(sn, 'nservers') else 1

        stage_spec['nodes'].append(node_spec)

    # Clone class specifications
    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
        njobs = sn.njobs[k] if hasattr(sn, 'njobs') else 1

        class_spec = {
            'id': k,
            'name': class_name,
            'type': 'open' if np.isinf(njobs) else 'closed',
            'population': njobs,
        }
        stage_spec['classes'].append(class_spec)

    # Create service processes with exponential rates
    for ist in range(sn.nstations):
        for k in range(sn.nclasses):
            sched = sn.sched[ist] if hasattr(sn, 'sched') else None
            sched_name = sched.name if hasattr(sched, 'name') else 'FCFS'

            # Skip source nodes
            if sched_name == 'EXT':
                # For source, keep original arrival rate
                if hasattr(sn, 'rates'):
                    rate = sn.rates[ist, k]
                    if not np.isnan(rate) and rate > 0:
                        stage_spec['processes'].append({
                            'station': ist,
                            'class': k,
                            'rate': rate,
                            'process_type': 'arrival',
                            'distribution': 'Exp',
                        })
            else:
                # For service, use the modulated exponential rate
                stage_spec['processes'].append({
                    'station': ist,
                    'class': k,
                    'rate': exp_rate,
                    'process_type': 'service',
                    'distribution': 'Exp',
                })

    return stage_spec


def _build_stage_network(model: Any, stage_name: str, exp_rate: float) -> Any:
    """
    Build a stage network by cloning the original model and replacing MMPP2 with Exp.

    Args:
        model: Original Network model
        stage_name: Name suffix for this stage network
        exp_rate: Exponential service rate to use (replacing MMPP)

    Returns:
        Network object for this stage
    """
    from line_solver import Network, Queue, Delay, Source, Sink
    from line_solver import OpenClass, ClosedClass
    from line_solver.distributions import Exp
    from line_solver.distributions.markovian import MMPP2

    # Create new network
    orig_name = model.name if hasattr(model, 'name') else 'model'
    stage_net = Network(f'{orig_name}_{stage_name}')

    # Map from original node names to new nodes
    node_map = {}

    # Get nodes from model
    orig_nodes = model.get_nodes() if hasattr(model, 'get_nodes') else (model.nodes if hasattr(model, 'nodes') else [])

    # PASS 1: Create all nodes
    for orig_node in orig_nodes:
        node_name = orig_node.name if hasattr(orig_node, 'name') else str(orig_node)
        node_class_name = orig_node.__class__.__name__

        if node_class_name == 'Source':
            new_node = Source(stage_net, node_name)
        elif node_class_name == 'Sink':
            new_node = Sink(stage_net, node_name)
        elif node_class_name == 'Delay':
            new_node = Delay(stage_net, node_name)
        elif node_class_name == 'Queue':
            from line_solver.lang import SchedStrategy
            sched = orig_node.sched_strategy if hasattr(orig_node, 'sched_strategy') else SchedStrategy.FCFS
            new_node = Queue(stage_net, node_name, sched)
            if hasattr(orig_node, 'num_servers') and orig_node.num_servers > 1:
                new_node.set_num_servers(orig_node.num_servers)
        else:
            continue

        node_map[node_name] = new_node

    # PASS 2: Create job classes
    orig_classes = model.get_classes() if hasattr(model, 'get_classes') else (model.classes if hasattr(model, 'classes') else [])
    for orig_class in orig_classes:
        class_name = orig_class.name if hasattr(orig_class, 'name') else str(orig_class)
        class_type = orig_class.__class__.__name__

        if class_type == 'OpenClass':
            OpenClass(stage_net, class_name)
        elif class_type == 'ClosedClass':
            # Get reference station - try method first, then attribute
            ref_stat = None
            if hasattr(orig_class, 'get_reference_station'):
                ref_stat = orig_class.get_reference_station()
            elif hasattr(orig_class, 'reference_station'):
                ref_stat = orig_class.reference_station

            # Get population - try method first, then attribute
            if hasattr(orig_class, 'get_population'):
                population = orig_class.get_population()
            elif hasattr(orig_class, 'population'):
                population = orig_class.population
            else:
                population = 1

            if ref_stat is not None:
                ref_name = ref_stat.name if hasattr(ref_stat, 'name') else str(ref_stat)
                if ref_name in node_map:
                    new_ref = node_map[ref_name]
                    ClosedClass(stage_net, class_name, int(population), new_ref)

    # Get stage network classes
    stage_classes = stage_net.get_classes() if hasattr(stage_net, 'get_classes') else []

    # PASS 3: Set arrivals from Source nodes
    for orig_node in orig_nodes:
        if orig_node.__class__.__name__ == 'Source' and orig_node.name in node_map:
            new_source = node_map[orig_node.name]
            for i, orig_class in enumerate(orig_classes):
                if i < len(stage_classes):
                    new_class = stage_classes[i]
                    if hasattr(orig_node, 'get_arrival'):
                        arr_dist = orig_node.get_arrival(orig_class)
                        if arr_dist is not None:
                            new_source.set_arrival(new_class, arr_dist)

    # PASS 4: Set services, replacing MMPP2 with Exp
    for orig_node in orig_nodes:
        node_class_name = orig_node.__class__.__name__
        if node_class_name in ('Queue', 'Delay') and orig_node.name in node_map:
            new_node = node_map[orig_node.name]
            for i, orig_class in enumerate(orig_classes):
                if i < len(stage_classes):
                    new_class = stage_classes[i]
                    if hasattr(orig_node, 'get_service'):
                        orig_dist = orig_node.get_service(orig_class)
                        if orig_dist is not None:
                            # Replace MMPP2 with Exp
                            if isinstance(orig_dist, MMPP2) or 'MMPP' in orig_dist.__class__.__name__:
                                new_dist = Exp(exp_rate)
                            else:
                                new_dist = orig_dist
                            new_node.set_service(new_class, new_dist)

    # PASS 5: Setup routing
    route_nodes = [node_map[n.name] for n in orig_nodes if n.name in node_map]
    if len(route_nodes) >= 2:
        try:
            stage_net.link(Network.serial_routing(*route_nodes))
        except Exception:
            pass

    return stage_net


__all__ = [
    'qn2line',
    'line2qn',
    'qn2lqn',
    'lqn2qn',
    'LQNTask',
    'LQNProcessor',
    'LQNModel',
    'mapqn2renv',
    'RandomEnvironmentModel',
    'MMPP2Params',
]
