"""
LINE Solver I/O Functions

This module provides I/O functions for exporting LINE network models to various formats.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from .indexed_table import IndexedTable
from .model_adapter import ModelAdapter

__all__ = ['IndexedTable', 'ModelAdapter', 'M2M', 'QN2JSIMG', 'qn2jsimg', 'LQN2QN', 'lqn2qn']


class M2M:
    """
    Model-to-Model conversion utility class.

    Provides methods for loading models from various formats (JSIM, MAT, etc.)
    and converting between model representations.
    """

    def JSIM2LINE(self, filename, modelName=None):
        """
        Load a JMT JSIM model file and convert to LINE Network.

        Args:
            filename: Path to the JSIM XML file
            modelName: Optional model name (default: extracted from file)

        Returns:
            Network: LINE network model

        Example:
            >>> m2m = M2M()
            >>> model = m2m.JSIM2LINE('model.jsimg')
        """
        from ..api.io.jmt_io import jsim2line
        from ..lang import Network, Queue, Delay, Source, Sink, Router
        from ..lang import OpenClass, ClosedClass
        from ..distributions import Exp, Erlang, HyperExp, Det
        from ..scheduling import SchedStrategy

        # Parse JSIM file
        spec = jsim2line(filename)

        # Create network
        name = modelName if modelName else spec.get('name', 'imported_model')
        model = Network(name)

        # Create nodes from specification
        node_map = {}
        for node_spec in spec.get('nodes', []):
            node_name = node_spec.get('name', f'Node{len(node_map)}')
            node_type = node_spec.get('type', 'Queue')

            if node_type == 'Source':
                node = Source(model, node_name)
            elif node_type == 'Sink':
                node = Sink(model, node_name)
            elif node_type == 'Delay':
                node = Delay(model, node_name)
            elif node_type == 'Router':
                node = Router(model, node_name)
            else:
                # Default to Queue with scheduling strategy
                sched = node_spec.get('scheduling', 'FCFS')
                sched_strategy = getattr(SchedStrategy, sched, SchedStrategy.FCFS)
                node = Queue(model, node_name, sched_strategy)
                if 'servers' in node_spec:
                    node.setNumberOfServers(node_spec['servers'])
                if 'capacity' in node_spec:
                    node.setCapacity(node_spec['capacity'])

            node_map[node_name] = node

        # Create classes from specification
        class_map = {}
        for class_spec in spec.get('classes', []):
            class_name = class_spec.get('name', f'Class{len(class_map)}')
            class_type = class_spec.get('type', 'closed')
            priority = class_spec.get('priority', 0)

            if class_type.lower() == 'open':
                job_class = OpenClass(model, class_name, priority)
            else:
                population = class_spec.get('population', 1)
                ref_station = class_spec.get('reference_station')
                ref_node = node_map.get(ref_station) if ref_station else None
                if ref_node is None:
                    # Use first delay or queue as reference
                    for n in model.nodes:
                        if isinstance(n, (Delay, Queue)):
                            ref_node = n
                            break
                job_class = ClosedClass(model, class_name, population, ref_node, priority)

            class_map[class_name] = job_class

        # Set service distributions
        for node_spec in spec.get('nodes', []):
            node_name = node_spec.get('name')
            node = node_map.get(node_name)
            if node is None:
                continue

            services = node_spec.get('services', {})
            for class_name, svc_spec in services.items():
                job_class = class_map.get(class_name)
                if job_class is None:
                    continue

                # Create distribution from service spec
                dist = self._create_distribution(svc_spec)
                if dist is not None:
                    if hasattr(node, 'setService'):
                        node.setService(job_class, dist)
                    elif hasattr(node, 'setArrival'):
                        node.setArrival(job_class, dist)

        # Set routing from connections
        routing = spec.get('routing', {})
        for (from_name, to_name), prob in routing.items():
            from_node = node_map.get(from_name)
            to_node = node_map.get(to_name)
            if from_node and to_node:
                model.addLink(from_node, to_node)

        return model

    def _create_distribution(self, svc_spec):
        """Create a distribution from a service specification."""
        from ..distributions import Exp, Erlang, HyperExp, Det, Disabled

        if svc_spec is None:
            return None

        dist_type = svc_spec.get('type', 'Exponential')
        params = svc_spec.get('params', {})

        if dist_type == 'Disabled':
            return Disabled.getInstance()
        elif dist_type == 'Exponential':
            rate = params.get('rate', params.get('lambda', 1.0))
            return Exp(rate)
        elif dist_type == 'Erlang':
            rate = params.get('rate', 1.0)
            k = params.get('k', params.get('shape', 1))
            return Erlang.fitMeanAndOrder(1.0 / rate, k)
        elif dist_type == 'Deterministic':
            value = params.get('value', params.get('t', 1.0))
            return Det(value)
        elif dist_type == 'Hyperexponential':
            p = params.get('p', 0.5)
            lambda1 = params.get('lambda1', 1.0)
            lambda2 = params.get('lambda2', 2.0)
            return HyperExp(p, lambda1, lambda2)
        else:
            # Default to exponential with rate 1
            return Exp(1.0)

    def MAT2LINE(self, filename):
        """
        Load a MATLAB .mat file and convert to LINE Network.

        Args:
            filename: Path to the .mat file

        Returns:
            Network: LINE network model

        Raises:
            NotImplementedError: MAT loading not yet implemented
        """
        raise NotImplementedError(
            "MAT2LINE requires scipy.io.loadmat and model parsing implementation. "
            "Currently not available in pure Python mode."
        )


def QN2JSIMG(model, outputFileName=None, options=None):
    """
    Writes a Network model to JMT JSIMG format.

    Creates a JSIMG (JMT simulation) XML file from a Network model.
    This file can be opened in JMT's graphical editor or used for simulation.

    Args:
        model: Network model to export
        outputFileName: Optional output file path (default: temp file)
        options: Optional solver options dictionary

    Returns:
        Path to the created JSIMG file

    Example:
        >>> model = Network('example')
        >>> # ... define model ...
        >>> fname = QN2JSIMG(model)
        >>> jsimgView(fname)  # View in JMT
    """
    from ..api.io.jmt_io import qn2jsimg as _qn2jsimg
    return _qn2jsimg(model, outputFileName, options)


# Alias for snake_case
qn2jsimg = QN2JSIMG


def LQN2QN(lqn):
    """
    Convert a LayeredNetwork (LQN) to a Network (QN) using REPLY signals.

    This converter transforms a Layered Queueing Network into an equivalent
    Queueing Network that uses REPLY signals to model synchronous call blocking.

    REPLY Signal Semantics:
    - Each synchCall creates a Request class and a Reply signal class
    - The caller queue blocks after completing service until Reply arrives
    - The callee processes the request and class-switches to Reply on completion
    - Reply signal unblocks the caller and continues downstream

    Network Topology:
      Think -> CallerQueue -> CalleeQueue -> CallerQueue (reply) -> Think
                    | blocks      | class switch to Reply    | unblocks

    This approach:
    - Correctly models BLOCKING semantics (server waits for reply)
    - Provides per-task queue metrics (queue length, utilization)
    - Supports multi-tier call chains (A -> B -> C)
    - Uses DES solver with REPLY signal support

    Args:
        lqn: LayeredNetwork model to convert.

    Returns:
        Network: Equivalent queueing network approximating the LQN behavior.

    Example:
        >>> lqn = LayeredNetwork('MyLQN')
        >>> # ... define LQN model ...
        >>> model = LQN2QN(lqn)
        >>> SolverSSA(model).getAvgTable()
    """
    from ..lang import Network, Queue, Delay, ClosedClass, Signal, SignalType
    from ..distributions import Exp
    from ..scheduling import SchedStrategy, RoutingStrategy
    from ..layered import CallType

    # Get LQN structure
    lsn = lqn.getStruct()

    # Create QN model
    model = Network(f"{lqn.name}-QN")

    # Identify reference tasks
    ref_task_indices = []
    for t in range(1, lsn.ntasks + 1):
        tidx = lsn.tshift + t
        if lsn.isref[tidx, 0] == 1:
            ref_task_indices.append(tidx)

    n_ref_tasks = len(ref_task_indices)
    if n_ref_tasks == 0:
        raise ValueError("LQN must have at least one reference task.")

    # Build task service demands
    task_service = {}
    for t in range(1, lsn.ntasks + 1):
        tidx = lsn.tshift + t
        # Sum up service demands from all activities of this task
        total_demand = 0.0
        if tidx in lsn.actsof:
            for aidx in lsn.actsof[tidx]:
                if aidx in lsn.hostdem:
                    demand = lsn.hostdem[aidx]
                    if demand and demand > 0:
                        total_demand += demand
        task_service[tidx] = total_demand

    # Build call graph: for each task, find all tasks it calls synchronously
    synch_calls_from = {}
    for t in range(1, lsn.ntasks + 1):
        tidx = lsn.tshift + t
        calls = []
        if tidx in lsn.actsof:
            for aidx in lsn.actsof[tidx]:
                # Check if this activity has calls
                if lsn.callpair is not None:
                    for c in range(1, lsn.ncalls + 1):
                        if lsn.callpair[c, 1] == aidx:
                            target_eidx = int(lsn.callpair[c, 2])
                            target_tidx = int(lsn.parent[target_eidx, 0])
                            call_mean = lsn.callpair[c, 3] if lsn.callpair.shape[1] > 3 else 1.0
                            # Check if sync call
                            if lsn.issynccaller[tidx, target_tidx] > 0:
                                calls.append({
                                    'target_tidx': target_tidx,
                                    'target_eidx': target_eidx,
                                    'call_mean': call_mean
                                })
        synch_calls_from[tidx] = calls

    # Find all tasks in call chains starting from reference tasks
    def collect_tasks_in_chain(start_tidx, visited=None):
        if visited is None:
            visited = set()
        if start_tidx in visited:
            return visited
        visited.add(start_tidx)
        for call in synch_calls_from.get(start_tidx, []):
            collect_tasks_in_chain(call['target_tidx'], visited)
        return visited

    tasks_in_chain = set()
    for ref_tidx in ref_task_indices:
        tasks_in_chain.update(collect_tasks_in_chain(ref_tidx))

    # Create nodes for all tasks
    task_queues = {}
    think_nodes = {}

    # Create think nodes for reference tasks
    for ref_tidx in ref_task_indices:
        task_name = lsn.names[ref_tidx]
        think_node = Delay(model, f"{task_name}_Think")
        think_nodes[ref_tidx] = think_node

    # Create queues for all tasks in call chains
    for tidx in tasks_in_chain:
        task_name = lsn.names[tidx]
        proc_idx = int(lsn.parent[tidx, 0]) if lsn.parent[tidx, 0] > 0 else 0
        n_servers = int(lsn.mult[0, proc_idx]) if proc_idx > 0 else 1
        sched = lsn.sched.get(proc_idx, SchedStrategy.FCFS)

        if n_servers >= 1e9 or sched == SchedStrategy.INF:
            queue = Delay(model, task_name)
        else:
            queue = Queue(model, task_name, sched if isinstance(sched, SchedStrategy) else SchedStrategy.FCFS)
            queue.setNumberOfServers(max(1, n_servers))
        task_queues[tidx] = queue

    # Create job classes for each reference task
    request_classes = {}
    reply_signals = {}

    for ref_tidx in ref_task_indices:
        task_name = lsn.names[ref_tidx]
        think_node = think_nodes[ref_tidx]
        population = int(lsn.mult[0, ref_tidx])

        # Create request class (closed)
        request_class = ClosedClass(model, f"{task_name}_Req", population, think_node)
        request_classes[ref_tidx] = request_class

        # Create reply signal class
        reply_signal = Signal(model, f"{task_name}_Reply", SignalType.REPLY)
        reply_signal.forJobClass(request_class)
        reply_signals[ref_tidx] = reply_signal

        # Override default RAND routing - set DISABLED so link() can set proper routes
        for node in model.nodes:
            if not hasattr(node, '__class__') or node.__class__.__name__ != 'Sink':
                node.setRouting(reply_signal, RoutingStrategy.DISABLED)

    # Set service times for all nodes and classes
    FINE_TOL = 1e-12
    for ref_tidx in ref_task_indices:
        request_class = request_classes[ref_tidx]
        reply_signal = reply_signals[ref_tidx]
        think_node = think_nodes[ref_tidx]

        # Think time at think node
        think_mean = lsn.think.get(ref_tidx, 0.0)
        if think_mean is None or think_mean < FINE_TOL:
            think_node.setService(request_class, Exp(1e8))
        else:
            think_node.setService(request_class, Exp(1.0 / think_mean))
        # Reply passes through think instantly
        think_node.setService(reply_signal, Exp(1e9))

        # Set service times at all task queues
        for tidx in tasks_in_chain:
            if tidx not in task_queues:
                continue
            queue = task_queues[tidx]
            service_mean = task_service.get(tidx, 0.0)

            if service_mean > FINE_TOL:
                queue.setService(request_class, Exp(1.0 / service_mean))
            else:
                queue.setService(request_class, Exp(1e8))
            # Reply signal passes through instantly
            queue.setService(reply_signal, Exp(1e9))

    # Build routing matrix with class switching for REPLY signals
    P = model.initRoutingMatrix()

    def build_call_chain(start_tidx):
        """Build the call chain from a starting task to the leaf."""
        chain = [start_tidx]
        current = start_tidx
        while True:
            calls = synch_calls_from.get(current, [])
            if not calls:
                break
            next_tidx = calls[0]['target_tidx']
            if next_tidx in chain:
                break  # Avoid cycles
            chain.append(next_tidx)
            current = next_tidx
        return chain

    for ref_tidx in ref_task_indices:
        request_class = request_classes[ref_tidx]
        reply_signal = reply_signals[ref_tidx]
        think_node = think_nodes[ref_tidx]

        # Build the complete call chain from reference task to leaf
        call_chain = build_call_chain(ref_tidx)

        if len(call_chain) == 1:
            # No synch calls - just loop Think -> RefQueue -> Think
            if ref_tidx in task_queues:
                ref_queue = task_queues[ref_tidx]
                P.set(request_class, request_class, think_node, ref_queue, 1.0)
                P.set(request_class, request_class, ref_queue, think_node, 1.0)
            else:
                P.set(request_class, request_class, think_node, think_node, 1.0)
            continue

        # Build Request path: Think -> task1 -> task2 -> ... -> leafTask
        first_tidx = call_chain[0]
        if first_tidx in task_queues:
            first_queue = task_queues[first_tidx]
            P.set(request_class, request_class, think_node, first_queue, 1.0)
        elif len(call_chain) > 1:
            first_queue = task_queues[call_chain[1]]
            P.set(request_class, request_class, think_node, first_queue, 1.0)
            call_chain = call_chain[1:]

        # Request path through the call chain
        for c in range(len(call_chain) - 1):
            from_tidx = call_chain[c]
            to_tidx = call_chain[c + 1]
            if from_tidx in task_queues and to_tidx in task_queues:
                from_queue = task_queues[from_tidx]
                to_queue = task_queues[to_tidx]
                P.set(request_class, request_class, from_queue, to_queue, 1.0)

        # Class switch at the LEAF node (last in chain)
        leaf_tidx = call_chain[-1]
        leaf_queue = task_queues[leaf_tidx]

        if len(call_chain) == 1:
            # Only one task in chain
            P.set(request_class, reply_signal, leaf_queue, think_node, 1.0)
        else:
            # Class switch at leaf, reply flows back through chain
            prev_tidx = call_chain[-2]
            prev_queue = task_queues[prev_tidx]
            P.set(request_class, reply_signal, leaf_queue, prev_queue, 1.0)

            # Reply flows back through intermediate nodes
            for c in range(len(call_chain) - 2, 0, -1):
                from_tidx = call_chain[c]
                to_tidx = call_chain[c - 1]
                from_queue = task_queues[from_tidx]
                to_queue = task_queues[to_tidx]
                P.set(reply_signal, reply_signal, from_queue, to_queue, 1.0)

            # Final hop: first task in chain -> Think (class switch back to Request)
            first_queue = task_queues[call_chain[0]]
            P.set(reply_signal, request_class, first_queue, think_node, 1.0)

    model.link(P)

    return model


# Alias for snake_case
lqn2qn = LQN2QN
