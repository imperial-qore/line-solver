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


def qn2lqn(model: Any) -> Any:
    """
    Convert a Queueing Network to Layered Queueing Network representation.

    Creates a LayeredNetwork model from a QN model by mapping:
    - A pseudo host with INF servers
    - Per-chain reference tasks with entries
    - Queue/Delay nodes to hosts/tasks/entries/activities
    - Routing matrix to OR-fork activity precedences

    Args:
        model: Network model (must have getStruct() method)

    Returns:
        LayeredNetwork object that can be solved by SolverLQNS

    References:
        MATLAB: matlab/src/io/QN2LQN.m
    """
    from ...layered import (
        LayeredNetwork, Processor, Task, Entry, Activity,
        ActivityPrecedence, PrecedenceType
    )
    from ...constants import SchedStrategy
    from ...distributions import Immediate

    sn = model.getStruct() if hasattr(model, 'getStruct') else model.get_struct()

    nNodes = sn.nnodes
    nClasses = sn.nclasses
    nChains = sn.nchains

    model_name = model.getName() if hasattr(model, 'getName') else (model.get_name() if hasattr(model, 'get_name') else 'model')

    lqn = LayeredNetwork(model_name)

    # Pseudo host with INF servers, INF scheduling
    PH = Processor(lqn, model_name, np.inf, SchedStrategy.INF)

    # Reference tasks and entries per chain
    RT = [None] * nChains
    RE = [None] * nChains
    for c in range(nChains):
        inchain = sn.inchain[c]
        if hasattr(inchain, 'flatten'):
            inchain = inchain.flatten().astype(int)
        total_jobs = sum(int(sn.njobs[r]) for r in inchain)
        RT[c] = Task(lqn, f'RefTask_{c+1}', total_jobs, SchedStrategy.REF)
        RT[c].on(PH)
        RE[c] = Entry(lqn, f'Chain_{c+1}')
        RE[c].on(RT[c])

    # Create Hosts, Tasks, Entries, Activities for Queue/Delay nodes
    P = [None] * nNodes
    T = [None] * nNodes
    E = [[None] * nClasses for _ in range(nNodes)]
    A = [[None] * nClasses for _ in range(nNodes)]

    from ..sn import NodeType as NT

    for i in range(nNodes):
        node_type = int(sn.nodetype[i])
        if node_type == NT.QUEUE or node_type == NT.DELAY:
            ist = int(sn.nodeToStation[i])
            nservers = int(sn.nservers[ist]) if not np.isinf(sn.nservers[ist]) else np.inf
            sched_strat = sn.sched[ist] if ist in sn.sched else SchedStrategy.FCFS
            P[i] = Processor(lqn, sn.nodenames[i], nservers, sched_strat)
            T[i] = Task(lqn, f'T_{sn.nodenames[i]}', np.inf, SchedStrategy.INF)
            T[i].on(P[i])

            for r in range(nClasses):
                c = _find_chain(sn, r)
                visits = sn.visits[c]
                if visits[i, r] > 0:
                    E[i][r] = Entry(lqn, f'E{i+1}_{r+1}')
                    E[i][r].on(T[i])
                    # Get service process from the model
                    nodes = model.get_nodes()
                    classes = model.get_classes()
                    service_dist = nodes[i].get_service(classes[r])
                    A[i][r] = Activity(lqn, f'Q{i+1}_{r+1}', service_dist)
                    A[i][r].on(T[i]).bound_to(E[i][r]).replies_to(E[i][r])
        elif node_type == NT.CLASSSWITCH:
            pass  # no-op
        elif node_type != NT.SOURCE and node_type != NT.SINK:
            pass  # Skip unsupported types silently for now

    # Create pseudo-activities on reference tasks
    PA = [[[None] * nClasses for _ in range(nNodes)] for _ in range(nChains)]
    boundToRE = [None] * nChains  # (node_idx, class_idx) or None

    for i in range(nNodes):
        node_type = int(sn.nodetype[i])
        if node_type == NT.CLASSSWITCH:
            for r in range(nClasses):
                c = _find_chain(sn, r)
                if _has_incoming_routing(sn, i, r):
                    PA[c][i][r] = Activity(lqn, f'CS_{c+1}_{i+1}_{r+1}', Immediate.getInstance())
                    PA[c][i][r].on(RT[c])
        elif node_type == NT.QUEUE or node_type == NT.DELAY:
            for r in range(nClasses):
                c = _find_chain(sn, r)
                visits = sn.visits[c]
                if visits[i, r] > 0:
                    inchain = sn.inchain[c]
                    if hasattr(inchain, 'flatten'):
                        inchain = inchain.flatten().astype(int)
                    first_class = int(inchain[0])
                    refstat = int(sn.refstat[first_class])
                    if i == refstat and r == first_class:
                        PA[c][i][r] = Activity(lqn, f'A{i+1}_{r+1}', Immediate.getInstance())
                        PA[c][i][r].on(RT[c]).bound_to(RE[c]).synch_call(E[i][r])
                        boundToRE[c] = (i, r)
                    else:
                        PA[c][i][r] = Activity(lqn, f'A{i+1}_{r+1}', Immediate.getInstance())
                        PA[c][i][r].on(RT[c]).synch_call(E[i][r])

    # Build OR-fork precedences from routing matrix
    for c in range(nChains):
        inchain = sn.inchain[c]
        if hasattr(inchain, 'flatten'):
            inchain = inchain.flatten().astype(int)

        for i in range(nNodes):
            node_type_i = int(sn.nodetype[i])
            if node_type_i in (NT.QUEUE, NT.DELAY, NT.CLASSSWITCH):
                for r in inchain:
                    r = int(r)
                    orfork_prec = []
                    orfork_prob = []

                    for j in range(nNodes):
                        node_type_j = int(sn.nodetype[j])
                        if node_type_j in (NT.QUEUE, NT.DELAY, NT.CLASSSWITCH):
                            for s in inchain:
                                s = int(s)
                                pr = sn.rtnodes[i * nClasses + r, j * nClasses + s]
                                if pr > 0 and _has_incoming_routing(sn, i, r):
                                    if boundToRE[c] is not None:
                                        if boundToRE[c][0] == j and boundToRE[c][1] == s:
                                            if PA[c][i][r] is not None:
                                                end_act = Activity(lqn,
                                                    f'End_{c+1}_{i+1}_{r+1}',
                                                    Immediate.getInstance())
                                                end_act.on(RT[c])
                                                orfork_prec.append(end_act)
                                                orfork_prob.append(pr)
                                        else:
                                            if PA[c][j][s] is not None:
                                                orfork_prec.append(PA[c][j][s])
                                                orfork_prob.append(pr)

                    if orfork_prec and PA[c][i][r] is not None:
                        RT[c].add_precedence(
                            ActivityPrecedence.OrFork(PA[c][i][r], orfork_prec, orfork_prob))

    return lqn


def _find_chain(sn, class_r: int) -> int:
    """Find chain index for a given class."""
    for c in range(sn.nchains):
        if sn.chains[c, class_r] > 0:
            return c
    return 0


def _has_incoming_routing(sn, node_i: int, class_r: int) -> bool:
    """Check if any routing leads to (node_i, class_r)."""
    col = node_i * sn.nclasses + class_r
    return np.any(sn.rtnodes[:, col] > 0)




def lqn2qn(lqn_model, replication='auto'):
    """
    Convert a LayeredNetwork into an equivalent Network by expanding each
    entry's activity subgraph into a step graph.

    This is the canonical converter, shared with ``line_solver.io.LQN2QN``;
    it is re-exported here so that the ``api.io`` namespace exposes the same
    implementation rather than a second, weaker one.

    Args:
        lqn_model: LayeredNetwork object
        replication: 'auto', 'materialize' or 'pool', as in ``io.LQN2QN``

    Returns:
        Network model

    References:
        MATLAB: matlab/src/io/LQN2QN.m
    """
    from ...io import LQN2QN
    return LQN2QN(lqn_model, replication=replication)


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
    """Random-environment image of a network with MAP/MMPP service or arrivals.

    Retained name for the transformation now implemented by :func:`map2renv`,
    which generalizes it from a single MMPP2 service process to any number of
    MAP, MMPP2 or MMAP arrival and service processes of arbitrary phase order.

    Args:
        model: Network with at least one MAP/MMPP2/MMAP process
        options: Optional solver options (config['map_env_maxstages'] caps the stages)

    Returns:
        Environment model with exponential rates modulated by the phases.
    """
    env, _ = map2renv(model, options)
    return env


def map2renv(model: Any, options: Optional[Dict] = None):
    """Markov-modulated image of a network as a queueing network in a random environment.

    Every non-renewal process is a point process modulated by the CTMC with
    generator Q = D0 + D1, whose conditional intensity in phase k is
    lambda(k) = sum_j D1(k,j). The transformation freezes each phase into an
    environment stage in which the process is the Poisson process of that
    intensity, i.e. an exponential arrival or service time, and lets the
    environment switch stages at the rates of Q. With P modulated processes the
    stage set is the Cartesian product of their phase spaces and the environment
    generator is the Kronecker sum of the individual Q's, so only one process
    changes phase at a time, as in the original model.

    The image is exact in structure for an MMPP (diagonal D1); for a general MAP
    the phase jumps that occur AT an event epoch (off-diagonal D1) are aggregated
    into Q, so the image matches the modulating chain and the conditional
    intensities but not the full inter-event autocorrelation.

    Populations are carried across stage switches unchanged (identity reset), as
    a phase switch moves no job.

    Mirrors matlab/src/io/map2renv.m.

    Args:
        model: Network with at least one MAP/MMPP2/MMAP process
        options: Optional solver options; config['map_env_maxstages'] caps the
            number of stages (default 64)

    Returns:
        Tuple (env, info) where info holds nstages, orders, is_mmpp, max_hold_time
        and the modulation records.
    """
    from ...environment import Environment
    from ...distributions import Exp
    from ..sn import sn_map_modulation

    zero = 1e-14
    max_stages = 64
    cfg = getattr(options, 'config', None) if options is not None else None
    if isinstance(options, dict):
        cfg = options.get('config', options)
    if isinstance(cfg, dict) and cfg.get('map_env_maxstages'):
        max_stages = int(cfg['map_env_maxstages'])

    sn = model.getStruct() if hasattr(model, 'getStruct') else model.get_struct()
    mods = sn_map_modulation(sn)
    if not mods:
        raise RuntimeError('The model declares no MAP, MMPP2 or MMAP process, '
                           'so it has no random-environment image.')

    orders = [int(m['order']) for m in mods]
    nstages = int(np.prod(orders))
    is_mmpp = all(m['is_mmpp'] for m in mods)
    if nstages > max_stages:
        raise RuntimeError(
            'The random-environment image of this model has %d stages (phase orders %s), above the '
            "options.config['map_env_maxstages'] cap of %d. Reduce the order of the modulating "
            'processes or raise the cap.' % (nstages, orders, max_stages))

    # Stage s enumerates the phase tuples in column-major order.
    P = len(mods)
    phase_of = np.zeros((nstages, P), dtype=int)
    for s in range(nstages):
        rem = s
        for p in range(P):
            phase_of[s, p] = rem % orders[p]
            rem //= orders[p]

    base_name = model.getName() if callable(getattr(model, 'getName', None)) else getattr(model, 'name', 'model')
    env = Environment(base_name + '_renv', nstages)
    names = []
    for s in range(nstages):
        name = 'Phase' + ''.join('_%d' % (phase_of[s, p] + 1) for p in range(P))
        names.append(name)
        env.add_stage(s, name, 'item', _build_map_stage(model, mods, phase_of[s], name, zero))

    # Kronecker sum of the phase generators: one process changes phase at a time.
    exit_rate = np.zeros(nstages)
    for s in range(nstages):
        for p in range(P):
            Qp = mods[p]['D0'] + sum(mods[p]['D1'])
            k = int(phase_of[s, p])
            stride = int(np.prod(orders[:p])) if p > 0 else 1
            for l in range(orders[p]):
                if l == k or Qp[k, l] <= zero:
                    continue
                t = s + (l - k) * stride
                env.add_transition(s, t, Exp(float(Qp[k, l])))
                exit_rate[s] += Qp[k, l]

    env.init()
    positive = exit_rate[exit_rate > 0]
    max_hold = float(np.max(1.0 / positive)) if positive.size else 0.0
    info = {'nstages': nstages, 'orders': orders, 'is_mmpp': bool(is_mmpp),
            'max_hold_time': max_hold, 'mods': mods}
    return env, info


def _build_map_stage(model: Any, mods, phases, stage_name: str, zero: float):
    """Copy of the base model with every modulated process frozen to its
    phase-conditional exponential rate."""
    from ...distributions import Exp

    stage_net = model.copy()
    base_name = model.getName() if callable(getattr(model, 'getName', None)) else getattr(model, 'name', 'model')
    if callable(getattr(stage_net, 'setName', None)):
        stage_net.setName('%s_%s' % (base_name, stage_name))
    else:
        stage_net.name = '%s_%s' % (base_name, stage_name)
    stations = stage_net.get_stations()
    classes = stage_net.get_classes()
    for p, mod in enumerate(mods):
        station = stations[mod['ist']]
        for c, r in enumerate(mod['classes']):
            rate = float(np.sum(mod['D1'][c][int(phases[p]), :]))
            if mod['arrival']:
                # A silent phase (zero intensity) is an ON/OFF source: keep it as
                # a rate rather than a disabled class, so that the class still
                # exists in every stage and the rate-averaged limit averages a zero.
                station.set_arrival(classes[r], Exp(max(rate, zero)))
            else:
                if rate <= zero:
                    raise RuntimeError(
                        'Phase %d of the service process of class %d at station %d has zero completion rate: '
                        'the station never empties while the environment sits in that stage, so the stage has '
                        'no steady state and the random-environment image is not defined. Model the stalled '
                        'server as a breakdown stage instead.' % (int(phases[p]) + 1, r + 1, mod['ist'] + 1))
                station.set_service(classes[r], Exp(rate))
    stage_net.refresh_struct()
    return stage_net



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
    'mapqn2renv',
    'map2renv',
    'RandomEnvironmentModel',
    'MMPP2Params',
]
