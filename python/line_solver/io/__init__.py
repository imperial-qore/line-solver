"""
LINE Solver I/O Functions

This module provides I/O functions for exporting LINE network models to various formats.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from .indexed_table import IndexedTable
from .model_adapter import ModelAdapter
from .linemodel_io import save_model, load_model

__all__ = ['IndexedTable', 'ModelAdapter', 'M2M', 'QN2JSIMG', 'qn2jsimg', 'LQN2QN', 'lqn2qn',
           'save_model', 'load_model']


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
        from ..lang import Place, Transition
        from ..lang import OpenClass, ClosedClass
        from ..lang.base import DropStrategy
        from ..distributions import Disabled
        from ..constants import SchedStrategy, RoutingStrategy, TimingStrategy, GlobalConstants

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
            elif node_type == 'Place':
                # SPN Place: token store. Capacity/drop-rule/marking are applied
                # after the classes exist (see the SPN configuration block below).
                node = Place(model, node_name)
            elif node_type == 'Transition':
                # SPN Transition: firing modes are added after the classes exist.
                node = Transition(model, node_name)
            else:
                # Default to Queue with scheduling strategy
                sched = node_spec.get('scheduling', 'FCFS')
                sched_strategy = getattr(SchedStrategy, sched, SchedStrategy.FCFS)
                node = Queue(model, node_name, sched_strategy)
                if 'servers' in node_spec:
                    node.setNumberOfServers(node_spec['servers'])
                if 'capacity' in node_spec:
                    # see _kb/12-interfaces-and-docs.md (JSIM2LINE builder: JMT sentinel conventions) for rationale
                    _cap = node_spec['capacity']
                    if _cap is not None and _cap > 0:
                        node.setCapacity(_cap)

            node_map[node_name] = node

        # Create classes from specification. Open classes reference the Source
        # (referenceSource); closed classes reference their reference station.
        class_map = {}
        for class_spec in spec.get('classes', []):
            class_name = class_spec.get('name', f'Class{len(class_map)}')
            class_type = class_spec.get('type', 'closed')
            priority = class_spec.get('priority', 0)

            if class_type.lower() == 'open':
                job_class = OpenClass(model, class_name, priority)
            else:
                population = class_spec.get('population', 1)
                ref_station = class_spec.get('refstation')
                ref_node = node_map.get(ref_station) if ref_station else None
                if ref_node is None:
                    # Use first delay or queue as reference
                    for n in model.getNodes():
                        if isinstance(n, (Delay, Queue)):
                            ref_node = n
                            break
                job_class = ClosedClass(model, class_name, population, ref_node, priority)

            class_map[class_name] = job_class

        # Ordered class list (creation order == model class-index order), used to
        # place the initial SPN marking into a per-class vector below.
        ordered_classes = list(class_map.values())

        # see _kb/12-interfaces-and-docs.md (JSIM2LINE builder: JMT sentinel conventions) for rationale
        drop_map = {
            'drop': DropStrategy.DROP,
            'BAS blocking': DropStrategy.BAS,
            'waiting queue': DropStrategy.WaitingQueue,
        }
        for node_spec in spec.get('nodes', []):
            node = node_map.get(node_spec.get('name'))
            if node is None:
                continue

            if isinstance(node, Place):
                total_cap = node_spec.get('total_capacity')
                # JMT encodes an unbounded place as -1; LINE keeps the default Inf.
                if total_cap is not None and total_cap != -1:
                    node.setCapacity(total_cap)
                for cls_name, cap in node_spec.get('class_capacities', {}).items():
                    jc = class_map.get(cls_name)
                    if jc is not None and cap is not None and cap != -1:
                        node.setClassCapacity(jc, cap)
                for cls_name, rule in node_spec.get('drop_rules', {}).items():
                    jc = class_map.get(cls_name)
                    if jc is not None:
                        node.setDropRule(jc, drop_map.get(rule, DropStrategy.DROP))

            elif isinstance(node, Transition):
                for mode_spec in node_spec.get('modes', []):
                    mode = node.addMode(mode_spec.get('name', 'Mode'))

                    servers = mode_spec.get('servers', -1)
                    # JMT encodes infinite servers as -1.
                    if servers is None or int(servers) < 0:
                        node.setNumberOfServers(mode, GlobalConstants.MaxInt)
                    else:
                        node.setNumberOfServers(mode, int(servers))

                    kind, distr = mode_spec.get('timing', ('TIMED', None))
                    if kind == 'IMMEDIATE':
                        node.setTimingStrategy(mode, TimingStrategy.IMMEDIATE)
                    else:
                        node.setTimingStrategy(mode, TimingStrategy.TIMED)
                        dist = self._create_distribution(distr)
                        if dist is not None:
                            node.setDistribution(mode, dist)

                    fp = mode_spec.get('firing_priority')
                    if fp is not None:
                        node.setFiringPriorities(mode, int(fp))
                    fw = mode_spec.get('firing_weight')
                    if fw is not None:
                        node.setFiringWeights(mode, float(fw))

                    # Enabling input arcs (tokens required from a place). JMT
                    # encodes an unbounded requirement as -1 -> Inf.
                    for st_name, entries in mode_spec.get('enabling', {}).items():
                        tgt = node_map.get(st_name)
                        if tgt is None:
                            continue
                        for cls_name, val in entries.items():
                            jc = class_map.get(cls_name)
                            if jc is not None:
                                node.setEnablingConditions(
                                    mode, jc, tgt,
                                    float('inf') if val == -1 else val)

                    # see _kb/12-interfaces-and-docs.md (JSIM2LINE builder: JMT sentinel conventions) for rationale
                    for st_name, entries in mode_spec.get('inhibiting', {}).items():
                        tgt = node_map.get(st_name)
                        if tgt is None:
                            continue
                        for cls_name, val in entries.items():
                            jc = class_map.get(cls_name)
                            if jc is None:
                                continue
                            if val == -1 or val == 0:
                                node.setInhibitingConditions(mode, jc, tgt, float('inf'))
                            else:
                                node.setInhibitingConditions(mode, jc, tgt, val)

                    # Firing output arcs (tokens produced into a place/sink).
                    for st_name, entries in mode_spec.get('firing', {}).items():
                        tgt = node_map.get(st_name)
                        if tgt is None:
                            continue
                        for cls_name, val in entries.items():
                            jc = class_map.get(cls_name)
                            if jc is not None:
                                node.setFiringOutcome(
                                    mode, jc, tgt,
                                    float('inf') if val == -1 else val)

        # see _kb/12-interfaces-and-docs.md (JSIM2LINE builder: JMT sentinel conventions) for rationale
        for node_spec in spec.get('nodes', []):
            node = node_map.get(node_spec.get('name'))
            if node is None:
                continue

            for class_name, arv_spec in node_spec.get('arrivals', {}).items():
                job_class = class_map.get(class_name)
                if job_class is None or not hasattr(node, 'setArrival'):
                    continue
                dist = self._create_distribution(arv_spec)
                node.setArrival(job_class, dist if dist is not None
                                else Disabled.getInstance())

            for class_name, svc_spec in node_spec.get('services', {}).items():
                job_class = class_map.get(class_name)
                if job_class is None or not hasattr(node, 'setService'):
                    continue
                dist = self._create_distribution(svc_spec)
                node.setService(job_class, dist if dist is not None
                                else Disabled.getInstance())

        # Wire topology from the JMT <connection> elements.
        for (from_name, to_name) in spec.get('connections', []):
            from_node = node_map.get(from_name)
            to_node = node_map.get(to_name)
            if from_node is not None and to_node is not None:
                model.addLink(from_node, to_node)

        # Apply per-class routing strategies parsed from the Router sections.
        strategy_map = {
            'Random': RoutingStrategy.RAND,
            'Probabilities': RoutingStrategy.PROB,
            'Round Robin': RoutingStrategy.RROBIN,
            'Join the Shortest Queue (JSQ)': RoutingStrategy.JSQ,
            'Disabled': RoutingStrategy.DISABLED,
        }
        for node_spec in spec.get('nodes', []):
            node = node_map.get(node_spec.get('name'))
            if node is None or isinstance(node, Sink) or not hasattr(node, 'setRouting'):
                continue
            for class_name, rspec in node_spec.get('routing', {}).items():
                job_class = class_map.get(class_name)
                if job_class is None:
                    continue
                strat = strategy_map.get(rspec.get('strategy'), RoutingStrategy.RAND)
                if strat == RoutingStrategy.PROB:
                    node.setRouting(job_class, RoutingStrategy.PROB)
                    for dest_name, prob in rspec.get('dests', {}).items():
                        dest = node_map.get(dest_name)
                        if dest is not None:
                            node.setProbRouting(job_class, dest, prob)
                else:
                    node.setRouting(job_class, strat)

        # see _kb/12-interfaces-and-docs.md (JSIM2LINE builder: JMT sentinel conventions) for rationale
        for st_name, per_class in spec.get('preload', {}).items():
            node = node_map.get(st_name)
            if not isinstance(node, Place):
                continue
            marking = [0.0] * len(ordered_classes)
            for cls_name, pop in per_class.items():
                jc = class_map.get(cls_name)
                if jc is not None and jc in ordered_classes:
                    marking[ordered_classes.index(jc)] = float(pop)
            node.setState(marking)

        return model

    def _create_distribution(self, svc_spec):
        """Create a distribution from a parsed JMT distribution spec.

        svc_spec is {'type': <jmt distr name>, 'params': {name: float}} as
        produced by jmt_io._parse_jmt_distr, or None for an empty (disabled)
        strategy. The parameter names are the JMT distrPar field names.
        """
        from ..distributions import (Exp, Erlang, HyperExp, Det, Coxian, Gamma,
                                      Pareto, Weibull, Lognormal, Uniform, Disabled)

        if svc_spec is None:
            return None

        dist_type = svc_spec.get('type', 'Exponential')
        params = svc_spec.get('params', {})

        if dist_type == 'Disabled':
            return Disabled.getInstance()
        elif dist_type == 'Exponential':
            # JMT ExponentialPar stores the rate under 'lambda'.
            rate = params.get('lambda', params.get('rate', 1.0))
            return Exp(rate)
        elif dist_type == 'Erlang':
            # JMT ErlangPar: 'alpha' per-phase rate, 'r' number of phases.
            alpha = params.get('alpha', 1.0)
            r = int(params.get('r', params.get('k', 1)))
            return Erlang(alpha, r)
        elif dist_type == 'Deterministic':
            return Det(params.get('t', params.get('value', 1.0)))
        elif dist_type == 'Hyperexponential':
            return HyperExp(params.get('p', 0.5),
                            params.get('lambda1', 1.0),
                            params.get('lambda2', 2.0))
        elif dist_type == 'Coxian':
            # JMT CoxianPar: 'lambda0','lambda1' phase rates, 'p0' branch prob.
            l0 = params.get('lambda0', 1.0)
            l1 = params.get('lambda1', 1.0)
            p0 = params.get('p0', 1.0)
            return Coxian([1.0 / l0, 1.0 / l1], [p0, 1.0])
        elif dist_type == 'Gamma':
            return Gamma(params.get('alpha', 1.0), params.get('beta', 1.0))
        elif dist_type == 'Pareto':
            return Pareto(params.get('alpha', 2.0), params.get('k', 1.0))
        elif dist_type == 'Weibull':
            return Weibull(params.get('alpha', 1.0), params.get('r', 1.0))
        elif dist_type == 'Lognormal':
            return Lognormal(params.get('mu', 0.0), params.get('sigma', 1.0))
        elif dist_type == 'Uniform':
            return Uniform(params.get('min', 0.0), params.get('max', 1.0))
        else:
            # Fall back to exponential with unit rate for unsupported JMT distrs.
            return Exp(1.0)

    def MAT2LINE(self, filename):
        """
        Load a MATLAB .mat file and convert to LINE Network.

        Loads a .mat file containing a saved LINE NetworkStruct (from
        model.getStruct() in MATLAB) and reconstructs the Network model.
        Supports common fields: nodenames, classnames, nodetype, nservers,
        sched, rates, routing, njobs, connmatrix.

        Args:
            filename: Path to the .mat file

        Returns:
            Network: LINE network model

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If .mat file does not contain a valid LINE struct
        """
        import os
        try:
            from scipy.io import loadmat
        except ImportError:
            raise ImportError("MAT2LINE requires scipy: pip install scipy")

        if not os.path.exists(filename):
            raise FileNotFoundError(f"MAT file not found: {filename}")

        mat = loadmat(filename, squeeze_me=True, struct_as_record=True)

        # Find the struct variable (skip MATLAB metadata keys)
        sn_data = None
        for key in mat:
            if not key.startswith('__'):
                val = mat[key]
                if hasattr(val, 'dtype') and val.dtype.names:
                    sn_data = val.flat[0] if val.ndim > 0 else val
                    break
        if sn_data is None:
            raise ValueError("No struct variable found in .mat file")

        from ..lang import Network, Queue, Delay, Source, Sink
        from ..lang import OpenClass, ClosedClass
        from ..distributions import Exp
        from ..constants import SchedStrategy
        from ..lang.base import NodeType
        import numpy as np

        def _get(field, default=None):
            try:
                return sn_data[field]
            except (ValueError, KeyError, IndexError):
                return default

        nstations = int(_get('nstations', 0))
        nnodes = int(_get('nnodes', nstations))
        nclasses = int(_get('nclasses', 0))

        nodenames = _get('nodenames', [])
        classnames = _get('classnames', [])
        nodetypes = _get('nodetype', [])
        nservers_arr = _get('nservers', [])
        sched_arr = _get('sched', [])
        rates = _get('rates', None)
        njobs = _get('njobs', [])
        connmatrix = _get('connmatrix', None)

        # Flatten MATLAB arrays
        if hasattr(nodenames, 'flat'):
            nodenames = [str(n).strip() for n in np.atleast_1d(nodenames)]
        if hasattr(classnames, 'flat'):
            classnames = [str(n).strip() for n in np.atleast_1d(classnames)]
        nodetypes = np.atleast_1d(nodetypes).flatten() if nodetypes is not None else []
        njobs = np.atleast_1d(njobs).flatten() if njobs is not None else []

        # Create Network
        model_name = os.path.splitext(os.path.basename(filename))[0]
        model = Network(model_name)

        # Create nodes
        nodes = []
        station_idx = 0
        for i in range(nnodes):
            name = nodenames[i] if i < len(nodenames) else f'Node{i}'
            nt = int(nodetypes[i]) if i < len(nodetypes) else 0

            if nt == NodeType.SOURCE.value:
                nodes.append(Source(model, name))
            elif nt == NodeType.SINK.value:
                nodes.append(Sink(model, name))
            elif nt == NodeType.DELAY.value:
                nodes.append(Delay(model, name))
            elif nt == NodeType.QUEUE.value:
                sched = SchedStrategy.FCFS
                if sched_arr is not None:
                    s_arr = np.atleast_1d(sched_arr).flatten()
                    if station_idx < len(s_arr):
                        sched_val = int(s_arr[station_idx])
                        for ss in SchedStrategy:
                            if ss.value == sched_val:
                                sched = ss
                                break
                q = Queue(model, name, sched)
                if nservers_arr is not None:
                    ns_arr = np.atleast_1d(nservers_arr).flatten()
                    if station_idx < len(ns_arr) and np.isfinite(ns_arr[station_idx]):
                        q.set_number_of_servers(int(ns_arr[station_idx]))
                nodes.append(q)
                station_idx += 1
            else:
                # Default to Queue
                nodes.append(Queue(model, name, SchedStrategy.FCFS))
                station_idx += 1

        # Create classes
        classes = []
        source_node = None
        delay_node = None
        for n in nodes:
            if isinstance(n, Source):
                source_node = n
            elif isinstance(n, Delay):
                delay_node = n

        for r in range(nclasses):
            cname = classnames[r] if r < len(classnames) else f'Class{r}'
            nj = float(njobs[r]) if r < len(njobs) else float('inf')
            if np.isinf(nj):
                # Open class
                classes.append(OpenClass(model, cname))
            else:
                # Closed class
                ref = delay_node if delay_node else (nodes[0] if nodes else None)
                classes.append(ClosedClass(model, cname, int(nj), ref))

        # Set service rates from rates matrix
        if rates is not None:
            rates = np.atleast_2d(rates)
            stations = model.get_stations()
            for ist in range(min(rates.shape[0], len(stations))):
                for r in range(min(rates.shape[1], len(classes))):
                    rate_val = float(rates[ist, r])
                    if rate_val > 0 and np.isfinite(rate_val):
                        stations[ist].set_service(classes[r], Exp(rate_val))

        # Set routing from connection matrix
        if connmatrix is not None:
            conn = np.atleast_2d(connmatrix)
            routing_nodes = []
            for i in range(min(conn.shape[0], len(nodes))):
                for j in range(min(conn.shape[1], len(nodes))):
                    if conn[i, j] > 0:
                        routing_nodes.append((nodes[i], nodes[j]))
            if routing_nodes:
                model.link(Network.serialRouting(*[n for pair in routing_nodes for n in pair]))

        return model


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

    Python port of ``matlab/src/io/LQN2QN.m``; the two must stay in step.

    Construction:

    - One station per host processor (scheduling and multiplicity taken from
      the processor). Tasks sharing a processor share the station, as in the
      LQN semantics where the processor is the contended resource.
    - One Delay per reference task, holding its think time.
    - One closed class per step of the expanded activity graph. A step is an
      activity, or one call stage of an activity that issues synchronous calls.
      Steps carry population 0 except the think class, which carries the
      reference task multiplicity.
    - A synchronous call site blocks its caller: the step class has a REPLY
      signal bound to it, the token proceeds to the callee, and the callee's
      replying activity class-switches to that signal, which returns to the
      caller station and unblocks it.
    - Call multiplicity mean m is unrolled into floor(m) mandatory call stages
      plus, if m is not integer, one further stage entered with probability
      m-floor(m).
    - AND precedences become Fork and Join nodes, with one Router per branch
      since a Fork cannot switch class per output link. A branch tail issuing a
      synchronous call is given a merge step, so that it reaches the Join in an
      ordinary class rather than as a REPLY signal, which carries no
      forked-task identity.
    - A CacheTask becomes a Cache node. The activity bound to an ItemEntry is
      the read step and sits on that node; its two CacheAccess successors
      become the hit and the miss class, and since the class switch is made by
      the Cache node itself, the routes leaving it are written in the successor
      class.
    - An asynchronous call is lowered to a non-blocking visit: the caller does
      not hold its server, but is serialised behind the callee, since a closed
      network cannot create the second token a concurrent send would need.
    - Entry forwarding splits the reply exits of the forwarding entry: with the
      forwarding probability the request is handed to the target entry, which
      replies to the original caller.
    - The multiplicity of a non-reference task is its thread pool: at most that
      many requests may be inside the task at once, including its nested
      callees. It is enforced by a finite capacity region with one linear
      admission constraint per task, and the calls of such a task do not hold
      the caller's server.

    - An entry with an open arrival process receives requests from a Source:
      they traverse the entry's subgraph in open classes and leave through a
      Sink where the entry would reply. Calls on an open chain do not hold the
      caller's server; thread pools of the traversed tasks are still enforced.

    - Phase-2 activities, the successors of a replying activity, run after the
      reply: the replying step's exit routes back to the caller, and each of
      its service completions spawns the continuation at the host station
      (sn.classspawn). The spawned token walks the phase-2 subgraph holding
      only the task's own thread and is destroyed at the chain end, through a
      NEGATIVE signal that always misses on a closed chain, or through the
      Sink on an open one. A boundary that ends on a call site is
      normalised through a merge step at the host station; a phase 2
      that opens with an AND-fork spawns into an immediate head that feeds
      the Fork; at an AND-join branch tail the spawned token inherits the
      fork identity of the trigger and stands in for it at the Join; at a
      cache read the reply is emitted by an immediate trigger step per
      hit/miss outcome, whose completion spawns the matching branch
      continuation.

    - An AND-join quorum k of n is applied to the Join node in the class that
      entered the Fork; k equal to the branch count is the default
      wait-for-all and is left alone. An activity think time becomes an extra
      step on a shared ActivityThink delay, in series with the host demand, so
      the task keeps its thread for it while its processor is released.

    Not yet represented: delayed-hit retrieval on the cache miss path, the
    thread pool of a task with an internal AND-fork, task and processor
    replication with its fan-out, and the setup and delay-off times of a
    function task. Each is reported through a warning.

    Args:
        lqn: LayeredNetwork model to convert.

    Returns:
        Network: queueing network modelling the LQN with REPLY signal blocking.

    Example:
        >>> lqn = LayeredNetwork('MyLQN')
        >>> # ... define LQN model ...
        >>> model = LQN2QN(lqn)
        >>> SolverLDES(model).getAvgTable()
    """
    import warnings
    import numpy as np

    from ..lang import (Network, Queue, Delay, Cache, Fork, Join, Router,
                        Source, Sink, ClosedClass, OpenClass, Signal, SignalType,
                        ClosedSignal)
    from ..distributions import Distribution, Exp, Immediate
    from ..constants import SchedStrategy, RoutingStrategy
    from ..lang.base import JoinStrategy, ReplacementStrategy

    FINE_TOL = 1e-8
    MAXCALLSTAGES = 20          # guard against unrolling a huge call multiplicity
    ID_POST_AND = 12
    ID_PRE_AND = 2
    CALL_SYNC, CALL_ASYNC, CALL_FWD = 1, 2, 3

    lsn = lqn.getStruct()
    model = Network("%s-QN" % lqn.name)

    ref_task_indices = [lsn.tshift + t for t in range(1, lsn.ntasks + 1)
                        if lsn.isref[lsn.tshift + t, 0] == 1]

    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    open_entries = []
    arrival_map = getattr(lsn, 'arrival', None) or {}
    for e in range(lsn.eshift + 1, lsn.eshift + lsn.nentries + 1):
        arv = arrival_map.get(e)
        if arv is not None and arv.getMean() > FINE_TOL and arv.getMean() != float('inf'):
            open_entries.append(e)

    if not ref_task_indices and not open_entries:
        raise ValueError("LQN must have at least one reference task or open arrival.")

    def sched_of(idx):
        s = lsn.sched.get(idx, SchedStrategy.FCFS)
        if isinstance(s, str):
            for cand in SchedStrategy:
                if cand.value == s:
                    return cand
            return SchedStrategy.FCFS
        return s

    def demand_of(aidx):
        """Host demand of an activity as a distribution, or None if negligible."""
        proc = getattr(lsn, 'hostdem_proc', {}) or {}
        if aidx in proc and proc[aidx] is not None:
            d = proc[aidx]
            if not isinstance(d, Immediate) and d.getMean() > FINE_TOL:
                return d
            return None
        m = lsn.hostdem.get(aidx, 0.0)
        if m is None or m <= FINE_TOL:
            return None
        return Exp.fitMean(m)

    # ---------------------------------------------------- unsupported features
    if lsn.calltype is not None:
        if any(int(lsn.calltype[c]) == CALL_ASYNC for c in range(1, lsn.ncalls + 1)):
            warnings.warn("LQN2QN: asynchronous calls are represented as non-blocking "
                          "visits: the caller releases its server but remains serialised "
                          "behind the callee.")
    if lsn.hasretrieval is not None and lsn.hasretrieval.any():
        warnings.warn("LQN2QN: delayed-hit retrieval on the cache miss path is not "
                      "represented.")
    if getattr(lsn, 'repl', None) is not None:
        repl = np.asarray(lsn.repl).ravel()
        over = [i for i in range(1, min(len(repl), lsn.nhosts + lsn.ntasks + 1)) if repl[i] > 1]
        if over:
            warnings.warn("LQN2QN: replication of %s is not represented: the replicas are "
                          "collapsed into a single station and their fan-out is ignored."
                          % lsn.names[over[0]])
    if getattr(lsn, 'isfunction', None) is not None and np.asarray(lsn.isfunction).any():
        warnings.warn("LQN2QN: setup and delay-off times of function tasks are not "
                      "represented: the task is converted as an ordinary always-on station.")
    # ------------------------------- tasks whose multiplicity is a thread pool
    # Such a task holds one thread per request from entry to reply, also across
    # its nested synchronous calls. Its calls do not hold the caller's server,
    # and a finite capacity region caps the jobs across its step classes at the
    # task multiplicity instead. Excluded, with a warning: a task with an
    # internal AND-fork (its forked siblings would double-count the thread that
    # spawned them), and a task called from inside an AND-fork branch (the
    # fork-join transformation retags branch flows into auxiliary classes,
    # which would silently bypass the task's admission row).
    def task_has_and_fork(tidx):
        if lsn.actposttype is None:
            return False
        for a in range(lsn.ashift + 1, lsn.ashift + lsn.nacts + 1):
            if int(lsn.parent[a, 0]) == tidx and a < len(lsn.actposttype) and \
                    int(lsn.actposttype[a]) == ID_POST_AND:
                return True
        return False

    fcr_task = set()
    for t in range(1, lsn.ntasks + 1):
        tidx = lsn.tshift + t
        if lsn.isref[tidx, 0] != 0:
            continue
        if lsn.iscache is not None and lsn.iscache[tidx, 0] != 0:
            continue
        mult = lsn.mult[0, tidx]
        if mult >= 1e9 or sched_of(tidx) == SchedStrategy.INF:
            continue
        if task_has_and_fork(tidx):
            warnings.warn("LQN2QN: multiplicity of task %s is not enforced: an "
                          "AND-fork inside a task cannot be capped by a finite "
                          "capacity region, whose job count would double-count "
                          "the forked siblings." % lsn.names[tidx])
            continue
        fcr_task.add(tidx)

    branch_acts = []
    if lsn.actposttype is not None:
        for a in range(lsn.ashift + 1, lsn.ashift + lsn.nacts + 1):
            if a >= len(lsn.actposttype) or int(lsn.actposttype[a]) != ID_POST_AND:
                continue
            frontier = [a]
            while frontier:
                cur = frontier.pop(0)
                if cur in branch_acts:
                    continue
                branch_acts.append(cur)
                if lsn.actpretype is not None and cur < len(lsn.actpretype) and \
                        int(lsn.actpretype[cur]) == ID_PRE_AND:
                    continue    # branch tail: do not traverse past the join
                for s in range(lsn.ashift + 1, lsn.ashift + lsn.nacts + 1):
                    if lsn.graph[cur, s] != 0 and int(lsn.parent[s, 0]) == int(lsn.parent[cur, 0]):
                        frontier.append(s)
    if branch_acts and fcr_task:
        front = []
        for a in branch_acts:
            for c in lsn.callsof.get(a, []):
                front.append(int(lsn.parent[int(lsn.callpair[c, 2]), 0]))
        shadow = set()
        while front:
            t_ = front.pop(0)
            if t_ in shadow:
                continue
            shadow.add(t_)
            for a in range(lsn.ashift + 1, lsn.ashift + lsn.nacts + 1):
                if int(lsn.parent[a, 0]) != t_:
                    continue
                for c in lsn.callsof.get(a, []):
                    front.append(int(lsn.parent[int(lsn.callpair[c, 2]), 0]))
        for t_ in sorted(shadow & fcr_task):
            fcr_task.discard(t_)
            warnings.warn("LQN2QN: multiplicity of task %s is not enforced: it is "
                          "called from inside an AND-fork branch, whose flows the "
                          "fork-join transformation retags outside the admission "
                          "constraint." % lsn.names[t_])

    # ------------------------------------------- stations, one per host processor
    host_station = {}
    host_is_delay = {}
    for h in range(1, lsn.nhosts + 1):
        nservers = lsn.mult[0, h]
        sched = sched_of(h)
        if nservers >= 1e9 or sched == SchedStrategy.INF:
            host_station[h] = Delay(model, lsn.names[h])
            host_is_delay[h] = True
        else:
            q = Queue(model, lsn.names[h], sched)
            q.setNumberOfServers(max(1, int(nservers)))
            host_station[h] = q
            host_is_delay[h] = False

    think_node = {}
    for ref_tidx in ref_task_indices:
        think_node[ref_tidx] = Delay(model, "%s_Think" % lsn.names[ref_tidx])

    # ------------------------------- pass 1: expand the activity graph into steps
    step_aidx, step_host, step_svc, step_name = [], [], [], []
    step_blocks, step_is_think, step_ref_task = [], [], []
    step_node, step_class_owner = [], []

    # flow rows: (from_step, to_step, prob, from_is_signal, in_target_class)
    flow = []
    # reply rows: (callee_exit_step, caller_step, callee_exit_is_signal, prob)
    reply = []
    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    spawn_pairs = []
    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    ph2_exits = []
    entry_stack = []
    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    thread_stack = []
    # Thread-pool tasks holding a thread while a job is at each step.
    step_tasks = []
    cache_node_of = {}
    cache_wiring = []
    # [join_step, join_aidx]; the quorum is applied once the fork class exists.
    join_quorum = []
    # Shared INF station carrying the activity think times: the task keeps its
    # thread across a think time but its host processor is released.
    act_think_node = []

    def add_step(aidx, hidx, svc, name, blocks, isthink, ref_tidx):
        step_aidx.append(aidx)
        step_host.append(hidx)
        step_svc.append(svc)
        step_name.append(name)
        step_blocks.append(blocks)
        step_is_think.append(isthink)
        step_ref_task.append(ref_tidx)
        step_node.append(None)
        sid = len(step_aidx) - 1
        step_class_owner.append(sid)
        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
        step_tasks.append(frozenset(t for t in thread_stack if t in fcr_task))
        return sid

    def add_aux_step(node_obj, owner_step, name, ref_tidx):
        """A step on a Fork, Join or Router node: no station, no class of its own."""
        sid = add_step(0, 0, None, name, False, False, ref_tidx)
        step_node[sid] = node_obj
        step_class_owner[sid] = step_class_owner[owner_step]
        return sid

    def add_route(from_port, to_step, prob):
        flow.append((from_port[0], to_step, prob, from_port[1], False))

    def add_cache_route(from_step, to_step):
        """Leaving a Cache node: the switch is made by the node, so the route is
        declared in the class of the target step."""
        flow.append((from_step, to_step, 1.0, False, True))

    def get_cache_node(tidx):
        if tidx in cache_node_of:
            return cache_node_of[tidx]
        caps = lsn.itemcap.get(tidx)
        cap = int(caps[0]) if caps is not None and len(caps) else 1
        cnode = Cache(model, "%s_Cache" % lsn.names[tidx], int(lsn.nitems[tidx, 0]),
                      cap, ReplacementStrategy(int(lsn.replacestrat[tidx, 0])))
        cache_node_of[tidx] = cnode
        return cnode

    def is_and_fork(succ):
        if len(succ) < 2 or lsn.actposttype is None:
            return False
        return all(int(lsn.actposttype[s]) == ID_POST_AND for s in succ)

    def is_and_join_pre(aidx):
        if lsn.actpretype is None:
            return False
        return 0 < aidx < len(lsn.actpretype) and int(lsn.actpretype[aidx]) == ID_PRE_AND

    def count_and_join_branches(join_aidx):
        """Branch tails feeding an AND-join, i.e. its PRE_AND predecessors."""
        return sum(1 for p in range(lsn.graph.shape[0])
                   if p != join_aidx and lsn.graph[p, join_aidx] != 0 and is_and_join_pre(p))

    def apply_join_quorum(join_node, join_class, join_aidx):
        """A join whose quorum equals its branch count already waits for all
        branches, the default JoinStrategy.STD, so only a genuine quorum k < n
        switches the node to JoinStrategy.PARTIAL."""
        if join_node is None or join_class is None or lsn.actquorum is None \
                or not (0 < join_aidx < len(lsn.actquorum)):
            return
        quorum = int(lsn.actquorum[join_aidx])
        nbranches = count_and_join_branches(join_aidx)
        if quorum < 1 or nbranches < 1 or quorum >= nbranches:
            return
        join_node.set_strategy(join_class, JoinStrategy.PARTIAL)
        join_node.set_required(join_class, quorum)

    def act_think_of(aidx):
        """Think time of an activity, or None when it has none. It is a delay in
        series with the host demand, held at the activity's own task (the thread
        is kept) but with the host processor released, as in lqns."""
        t = (getattr(lsn, 'actthink', None) or {}).get(aidx)
        if t is None:
            return None
        if isinstance(t, Distribution):
            return t if not isinstance(t, Immediate) and t.getMean() > FINE_TOL else None
        return Exp.fitMean(float(t)) if float(t) > FINE_TOL else None

    def act_think_station():
        """Single INF station shared by every activity think time."""
        if not act_think_node:
            act_think_node.append(Delay(model, "ActivityThink"))
        return act_think_node[0]

    def call_mean(cidx):
        return float(lsn.callpair[cidx, 3]) if lsn.callpair.shape[1] > 3 else 1.0

    def call_stages(aidx):
        """Unrolls the synchronous and asynchronous calls of an activity."""
        stages = []
        for cidx in lsn.callsof.get(aidx, []):
            ctype = int(lsn.calltype[cidx])
            if ctype not in (CALL_SYNC, CALL_ASYNC):
                continue
            isasync = (ctype == CALL_ASYNC)
            target_eidx = int(lsn.callpair[cidx, 2])
            m = call_mean(cidx)
            nfull = int(m + FINE_TOL)
            frac = m - nfull
            if nfull > MAXCALLSTAGES:
                warnings.warn("LQN2QN: call multiplicity %g truncated to %d stages."
                              % (m, MAXCALLSTAGES))
                nfull, frac = MAXCALLSTAGES, 0.0
            stages.extend([(target_eidx, 1.0, isasync)] * nfull)
            if frac > FINE_TOL:
                stages.append((target_eidx, frac, isasync))
        return stages

    def forwarding_of(eidx):
        """Rows (target_eidx, probability) of the forwarding calls of an entry."""
        fwd = []
        if lsn.calltype is None:
            return fwd
        for cidx in range(1, lsn.ncalls + 1):
            if int(lsn.calltype[cidx]) != CALL_FWD or int(lsn.callpair[cidx, 1]) != eidx:
                continue
            p = min(max(call_mean(cidx), 0.0), 1.0)
            if p > FINE_TOL:
                fwd.append((int(lsn.callpair[cidx, 2]), p))
        return fwd

    def expand_entry(eidx, ref_tidx):
        """Expands the activity subgraph bound to an entry, in the current context.

        Returns (first_step, reply_exits, terminals), the exits being
        (step, is_signal, prob) triples.
        """
        if eidx in entry_stack:
            warnings.warn("LQN2QN: recursive call cycle at entry %s truncated."
                          % lsn.names[eidx])
            return None, [], []
        entry_stack.append(eidx)
        thread_stack.append(int(lsn.parent[eidx, 0]))
        try:
            local_acts = lsn.actsof.get(eidx, [])
            bound = None
            for a in local_acts:
                if lsn.graph[eidx, a] != 0:
                    bound = a
                    break
            if bound is None:
                return None, [], []

            first_step, reply_exits, terminals = _walk_entry(bound, eidx, ref_tidx, local_acts)

            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            fwd = forwarding_of(eidx)
            if fwd and reply_exits:
                own_ports = list(reply_exits)
                fwd_exits = []
                pforw = 0.0
                # The forwarder's thread is released at the handoff, so the
                # forwarded chain is expanded without it on the thread stack.
                fwd_thread = thread_stack.pop()
                for target_eidx, p in fwd:
                    f_first, f_replies, f_terms = expand_entry(target_eidx, ref_tidx)
                    if f_first is None:
                        continue
                    for op in own_ports:
                        add_route((op[0], op[1]), f_first, op[2] * p)
                    pforw += p
                    fwd_exits.extend(f_replies)
                    fwd_exits.extend(f_terms)
                thread_stack.append(fwd_thread)
                # What is left of each of this entry's own ports still replies.
                residual = max(0.0, 1.0 - pforw)
                reply_exits = [(s, sig, pr * residual) for (s, sig, pr) in reply_exits]
                reply_exits.extend(fwd_exits)
            return first_step, reply_exits, terminals
        finally:
            entry_stack.pop()
            thread_stack.pop()

    def _walk_entry(a0, eidx, ref_tidx, local_acts):
        """Walks the intra-task activity graph of one entry."""
        visited_entry = {}
        visited_exit = {}
        join_of = {}
        fork_owner_stack = []
        state = {'saw_reply': False}
        reply_exits, terminals = [], []

        def replies_here(aidx):
            if lsn.replygraph is None:
                return False
            a = aidx - lsn.ashift - 1
            e = eidx - lsn.eshift - 1
            if a < 0 or a >= lsn.replygraph.shape[0] or e < 0 or e >= lsn.replygraph.shape[1]:
                return False
            return lsn.replygraph[a, e] != 0

        def make_activity_steps(aidx, tidx, cache_node):
            """One step for the host demand, plus one per unrolled call stage."""
            hidx = int(lsn.parent[tidx, 0])

            if cache_node is not None:
                # A read step holds no demand and issues no call: the lookup is
                # instantaneous, the work is done on the hit or miss branch.
                read_step = add_step(aidx, hidx, None, lsn.names[aidx], False, False, ref_tidx)
                step_node[read_step] = cache_node
                if lsn.callsof.get(aidx):
                    warnings.warn("LQN2QN: calls issued by cache read activity %s are "
                                  "ignored." % lsn.names[aidx])
                if demand_of(aidx) is not None:
                    warnings.warn("LQN2QN: host demand of cache read activity %s is "
                                  "ignored." % lsn.names[aidx])
                visited_exit[aidx] = (read_step, False)
                return read_step

            svc = demand_of(aidx)
            stages = call_stages(aidx)
            entry_step = add_step(aidx, hidx, svc, lsn.names[aidx], False, False, ref_tidx)
            # cur is the port through which the activity is currently left. A
            # blocking call site is left through its reply signal.
            cur = (entry_step, False)

            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            think_dist = act_think_of(aidx)
            if think_dist is not None:
                think_step = add_step(aidx, hidx, think_dist, "%s_think" % lsn.names[aidx],
                                      False, False, ref_tidx)
                step_node[think_step] = act_think_station()
                add_route(cur, think_step, 1.0)
                cur = (think_step, False)

            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            tmult = float(lsn.mult[0, tidx])
            host_blocks = (not host_is_delay.get(hidx, False) and tidx not in fcr_task
                           and ref_tidx != 0 and tmult < 1e9
                           and sched_of(tidx) != SchedStrategy.INF)

            for k, (target_eidx, prob, isasync) in enumerate(stages):
                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                blocks = host_blocks and not isasync
                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                if isasync:
                    async_thread = thread_stack.pop()
                callee_first, callee_replies, callee_terms = expand_entry(target_eidx, ref_tidx)
                if isasync:
                    thread_stack.append(async_thread)
                if callee_first is None:
                    continue        # callee not expandable: drop the call, never block
                # A callee path that neither replies nor continues still holds a
                # token, so it returns to the caller like a reply would.
                callee_replies = list(callee_replies) + list(callee_terms)

                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                needs_merge = (prob < 1.0) or (k < len(stages) - 1) or \
                    (not blocks) or is_and_join_pre(aidx)
                nxt = None
                if needs_merge:
                    nxt = add_step(aidx, hidx, None, "%s_c%d_ret" % (lsn.names[aidx], k + 1),
                                   False, False, ref_tidx)
                    if not blocks:
                        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                        step_node[nxt] = Router(model, "%s_c%d_ret" % (lsn.names[aidx], k + 1))

                if blocks:
                    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                    if k == 0 and prob >= 1.0 and cur == (entry_step, False):
                        blk = entry_step
                        step_blocks[blk] = True
                    else:
                        blk = add_step(aidx, hidx, None, "%s_c%d" % (lsn.names[aidx], k + 1),
                                       True, False, ref_tidx)
                        add_route(cur, blk, prob)
                        if prob < 1.0:
                            add_route(cur, nxt, 1.0 - prob)
                    add_route((blk, False), callee_first, 1.0)
                    for (s, sig, pr) in callee_replies:
                        reply.append((s, blk, sig, pr))
                    if needs_merge:
                        add_route((blk, True), nxt, 1.0)
                        cur = (nxt, False)
                    else:
                        cur = (blk, True)
                else:
                    add_route(cur, callee_first, prob)
                    if prob < 1.0:
                        add_route(cur, nxt, 1.0 - prob)
                    for (s, sig, pr) in callee_replies:
                        add_route((s, sig), nxt, pr)
                    cur = (nxt, False)

            visited_exit[aidx] = cur
            return entry_step

        def walk(aidx):
            if aidx in visited_entry:
                return visited_entry[aidx]

            tidx = int(lsn.parent[aidx, 0])

            # The activity bound to an ItemEntry of a CacheTask is the read
            # step: it sits on the Cache node rather than on the processor.
            cache_node = None
            if lsn.iscache is not None and lsn.iscache[tidx, 0] != 0 and lsn.graph[eidx, aidx] != 0:
                cache_node = get_cache_node(tidx)

            entry_step = make_activity_steps(aidx, tidx, cache_node)
            exit_port = visited_exit[aidx]
            visited_entry[aidx] = entry_step

            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            if replies_here(aidx):
                state['saw_reply'] = True

            succ = [s for s in local_acts if lsn.graph[aidx, s] != 0]
            if not succ:
                terminals.append((exit_port[0], exit_port[1], 1.0))
                return entry_step

            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            if replies_here(aidx):
                if cache_node is not None and len(succ) >= 2:
                    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                    trig_h = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                      "%s_ph2h" % lsn.names[aidx], False, False, ref_tidx)
                    trig_m = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                      "%s_ph2m" % lsn.names[aidx], False, False, ref_tidx)
                    add_cache_route(entry_step, trig_h)
                    add_cache_route(entry_step, trig_m)
                    cache_wiring.append((cache_node, entry_step, trig_h, trig_m,
                                         lsn.itemproc.get(eidx), int(lsn.nitems[eidx, 0])))
                    reply_exits.append((trig_h, False, 1.0))
                    reply_exits.append((trig_m, False, 1.0))
                    saved_stack = list(thread_stack)
                    del thread_stack[:]
                    thread_stack.append(tidx)
                    n_t0 = len(terminals)
                    for hm in range(2):
                        s_entry = walk(succ[hm])
                        if step_node[s_entry] is not None:
                            hm_head = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                               "%s_ph2b%d" % (lsn.names[aidx], hm + 1),
                                               False, False, ref_tidx)
                            add_route((hm_head, False), s_entry, 1.0)
                            s_entry = hm_head
                        spawn_pairs.append((trig_h if hm == 0 else trig_m, s_entry))
                    for (ts, tsig, tpr) in terminals[n_t0:]:
                        ph2_exits.append((ts, tsig, tpr, ref_tidx))
                    del terminals[n_t0:]
                    del thread_stack[:]
                    thread_stack.extend(saved_stack)
                    return entry_step
                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                ok_ctx = (cache_node is None
                          and (not fork_owner_stack or is_and_join_pre(aidx)))
                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                if ok_ctx and (exit_port[1]
                               or isinstance(step_node[exit_port[0]], Router)):
                    trig = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                    "%s_ph2t" % lsn.names[aidx], False, False, ref_tidx)
                    add_route(exit_port, trig, 1.0)
                    exit_port = (trig, False)
                ph2_spawn = (ok_ctx and not exit_port[1]
                             and step_node[exit_port[0]] is None)
                if not ph2_spawn:
                    warnings.warn("LQN2QN: phase-2 activities of %s run before the "
                                  "reply: the boundary is not a station departure, "
                                  "a degenerate cache read, or mid-branch inside "
                                  "an AND-fork." % lsn.names[aidx])
                else:
                    reply_exits.append((exit_port[0], exit_port[1], 1.0))
                    saved_stack = list(thread_stack)
                    del thread_stack[:]
                    thread_stack.append(tidx)
                    n_t0 = len(terminals)
                    pos_succ = [s for s in succ if lsn.graph[aidx, s] > 0]
                    target = None
                    if is_and_fork(succ):
                        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                        head = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                        "%s_ph2" % lsn.names[aidx], False, False, ref_tidx)
                        wire_and_fork((head, False), succ, aidx)
                        target = head
                    elif is_and_join_pre(aidx):
                        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                        head = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                        "%s_ph2" % lsn.names[aidx], False, False, ref_tidx)
                        wire_and_join((head, False), succ[0])
                        target = head
                    elif len(pos_succ) == 1:
                        s_entry = walk(pos_succ[0])
                        if step_node[s_entry] is None:
                            target = s_entry
                    if target is None:
                        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                        head = add_step(aidx, int(lsn.parent[tidx, 0]), None,
                                        "%s_ph2" % lsn.names[aidx], False, False, ref_tidx)
                        for s2 in pos_succ:
                            add_route((head, False), walk(s2), lsn.graph[aidx, s2])
                        target = head
                    spawn_pairs.append((exit_port[0], target))
                    for (ts, tsig, tpr) in terminals[n_t0:]:
                        ph2_exits.append((ts, tsig, tpr, ref_tidx))
                    del terminals[n_t0:]
                    del thread_stack[:]
                    thread_stack.extend(saved_stack)
                    return entry_step

            if cache_node is not None:
                # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
                if len(succ) < 2:
                    warnings.warn("LQN2QN: cache read %s has no hit/miss pair; treated "
                                  "as an ordinary activity." % lsn.names[aidx])
                else:
                    h_entry = walk(succ[0])
                    m_entry = walk(succ[1])
                    add_cache_route(entry_step, h_entry)
                    add_cache_route(entry_step, m_entry)
                    cache_wiring.append((cache_node, entry_step, h_entry, m_entry,
                                         lsn.itemproc.get(eidx), int(lsn.nitems[eidx, 0])))
                    return entry_step

            if is_and_fork(succ):
                wire_and_fork(exit_port, succ, aidx)
                return entry_step

            if is_and_join_pre(aidx):
                wire_and_join(exit_port, succ[0])
                return entry_step

            for s in succ:
                p = lsn.graph[aidx, s]
                if p <= 0:
                    continue
                add_route(exit_port, walk(s), p)
            return entry_step

        def wire_and_fork(from_port, fsucc, aidx):
            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            fork_node = Fork(model, "Fork_%s" % lsn.names[aidx])
            fork_step = add_aux_step(fork_node, from_port[0],
                                     "Fork_%s" % lsn.names[aidx], ref_tidx)
            add_route(from_port, fork_step, 1.0)
            fork_owner_stack.append(fork_step)
            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            ordered = ([s for s in fsucc if branch_replies(s)]
                       + [s for s in fsucc if not branch_replies(s)])
            for b, s in enumerate(ordered):
                rname = "Fork_%s_%d" % (lsn.names[aidx], b + 1)
                router_step = add_aux_step(Router(model, rname), fork_step, rname, ref_tidx)
                add_route((fork_step, False), router_step, 1.0)
                add_route((router_step, False), walk(s), 1.0)
            fork_owner_stack.pop()

        def wire_and_join(from_port, join_aidx):
            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            if join_aidx in join_of:
                add_route(from_port, join_of[join_aidx], 1.0)
                return
            if not fork_owner_stack:
                warnings.warn("LQN2QN: AND-join at %s has no enclosing AND-fork; "
                              "branches are serialised." % lsn.names[join_aidx])
                add_route(from_port, walk(join_aidx), 1.0)
                return
            fork_owner = fork_owner_stack[-1]
            join_node = Join(model, "Join_%s" % lsn.names[join_aidx],
                             step_node[fork_owner])
            join_step = add_aux_step(join_node, fork_owner,
                                     "Join_%s" % lsn.names[join_aidx], ref_tidx)
            add_route(from_port, join_step, 1.0)
            add_route((join_step, False), walk(join_aidx), 1.0)
            join_of[join_aidx] = join_step
            join_quorum.append((join_step, join_aidx))

        def branch_replies(a0):
            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            stack = [a0]
            seen = set()
            while stack:
                a = stack.pop()
                if a in seen:
                    continue
                seen.add(a)
                if replies_here(a):
                    return True
                if is_and_join_pre(a):
                    continue
                stack.extend(s for s in local_acts if lsn.graph[a, s] != 0)
            return False

        first_step = walk(a0)

        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
        if state['saw_reply']:
            reply_exits = reply_exits + terminals
            terminals = []
        return first_step, reply_exits, terminals

    for ref_tidx in ref_task_indices:
        think_step = add_step(0, 0, None, "%s_Think" % lsn.names[ref_tidx],
                              False, True, ref_tidx)
        for eidx in lsn.entriesof.get(ref_tidx, []):
            first_step, reply_exits, terminals = expand_entry(eidx, ref_tidx)
            if first_step is None:
                continue
            add_route((think_step, False), first_step, 1.0)
            # A reference task has no caller: its replies and its dead ends both
            # close the cycle at the think delay.
            for (s, sig, pr) in list(reply_exits) + list(terminals):
                add_route((s, sig), think_step, pr)

    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    open_wiring = []
    src_node = snk_node = None
    if open_entries:
        src_node = Source(model, 'Source')
        snk_node = Sink(model, 'Sink')
    for eidx in open_entries:
        first_step, reply_exits, terminals = expand_entry(eidx, 0)
        if first_step is None:
            warnings.warn("LQN2QN: open arrival entry %s has no bound activity; "
                          "ignored." % lsn.names[eidx])
            continue
        open_wiring.append((eidx, first_step, list(reply_exits) + list(terminals)))

    # ------------------------------- pass 2: create classes and reply signals
    nsteps = len(step_aidx)
    step_class = [None] * nsteps
    step_signal = [None] * nsteps

    for i in range(nsteps):
        if step_class_owner[i] != i:
            # Fork, Join and Router steps carry the job through unchanged.
            continue
        ref_tidx = step_ref_task[i]
        if ref_tidx == 0:
            # A step of an open arrival chain travels in an open class.
            step_class[i] = OpenClass(model, step_name[i])
        else:
            population = int(lsn.mult[0, ref_tidx]) if step_is_think[i] else 0
            step_class[i] = ClosedClass(model, step_name[i], population, think_node[ref_tidx])
    for i in range(nsteps):
        step_class[i] = step_class[step_class_owner[i]]

    # AND-join quorum, in the class the siblings are matched in.
    for (js, ja) in join_quorum:
        apply_join_quorum(step_node[js], step_class[js], ja)

    for i in range(nsteps):
        if step_blocks[i]:
            sig = Signal(model, "%s_Reply" % step_name[i], SignalType.REPLY)
            sig.forJobClass(step_class[i])
            step_signal[i] = sig

    # A Signal installs RAND routing at every node; clear it so that link()
    # only honours the routes set below.
    for i in range(nsteps):
        if step_signal[i] is None:
            continue
        for node in model.getNodes():
            if node.__class__.__name__ != 'Sink':
                node.setRouting(step_signal[i], RoutingStrategy.DISABLED)

    # Spawn bindings for phase-2 continuations: each completion of the trigger
    # class injects a job of the target class at the same station.
    for (trig, tgt) in spawn_pairs:
        step_class[trig].set_spawn_class(step_class[tgt])

    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    ph2_dump_node = None
    ph2_destructor_of = {}
    for (_, _, _, rft) in ph2_exits:
        if rft <= 0 or rft in ph2_destructor_of:
            continue
        if ph2_dump_node is None:
            ph2_dump_node = Queue(model, "Ph2Sink", SchedStrategy.FCFS)
        sig = ClosedSignal(model, "Ph2End_%s" % lsn.names[rft],
                           SignalType.NEGATIVE, think_node[rft])
        for node in model.getNodes():
            if node.__class__.__name__ != 'Sink':
                node.setRouting(sig, RoutingStrategy.DISABLED)
        ph2_dump_node.setService(sig, Immediate())
        ph2_destructor_of[rft] = sig

    # ---------------------------------------------------- pass 3: service times
    def station_of(i):
        if step_node[i] is not None:
            return step_node[i]
        if step_is_think[i]:
            return think_node[step_ref_task[i]]
        return host_station[step_host[i]]

    for i in range(nsteps):
        if step_node[i] is not None:
            if step_class_owner[i] == i and act_think_node and step_node[i] is act_think_node[0]:
                act_think_node[0].setService(step_class[i], step_svc[i])
                continue
            # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
            if step_class_owner[i] == i and isinstance(step_node[i], Router):
                if step_ref_task[i] == 0:
                    # Open chain: no think delay exists, declare the pair at the
                    # caller's host station, which the class never visits.
                    host_station[step_host[i]].setService(step_class[i], Immediate())
                else:
                    think_node[step_ref_task[i]].setService(step_class[i], Immediate())
            continue
        if step_is_think[i]:
            ref_tidx = step_ref_task[i]
            tmean = lsn.think.get(ref_tidx, 0.0)
            tnode = think_node[ref_tidx]
            if tmean is None or tmean <= FINE_TOL:
                tnode.setService(step_class[i], Immediate())
            else:
                tnode.setService(step_class[i], Exp.fitMean(tmean))
        else:
            station = host_station[step_host[i]]
            station.setService(step_class[i], step_svc[i] if step_svc[i] is not None else Immediate())

    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    for i in range(nsteps):
        if step_signal[i] is None:
            continue
        for h in range(1, lsn.nhosts + 1):
            host_station[h].setService(step_signal[i], Immediate())
        for tn in think_node.values():
            tn.setService(step_signal[i], Immediate())

    # Cache read/hit/miss wiring, now that the classes exist.
    for (cnode, read_step, hit_step, miss_step, itemproc, nitems) in cache_wiring:
        # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
        cnode.setRead(step_class[read_step], itemproc)
        cnode.setHitClass(step_class[read_step], step_class[hit_step])
        cnode.setMissClass(step_class[read_step], step_class[miss_step])

    # --------------------------------------------------------- pass 4: routing
    P = model.initRoutingMatrix()
    for (i, j, p, from_is_signal, in_target_class) in flow:
        if in_target_class:
            P.set(step_class[j], step_class[j], station_of(i), station_of(j), p)
        elif from_is_signal:
            P.set(step_signal[i], step_class[j], station_of(i), station_of(j), p)
        else:
            P.set(step_class[i], step_class[j], station_of(i), station_of(j), p)
    for (i, owner, via_signal, p) in reply:
        # A nested call returns through its own reply signal, which
        # class-switches into the reply signal of the outer call site.
        src = step_signal[i] if via_signal else step_class[i]
        P.set(src, step_signal[owner], station_of(i), station_of(owner), p)

    # -------------------------- phase-2 chain ends: destroy the spawned token
    for (s, sig, pr, rft) in ph2_exits:
        if rft == 0:
            # Open chain: the spawned token leaves through the Sink.
            ecls = step_class[s]
            P.set(ecls, ecls, station_of(s), snk_node, pr)
        elif sig:
            P.set(step_signal[s], ph2_destructor_of[rft], station_of(s), ph2_dump_node, pr)
        else:
            P.set(step_class[s], ph2_destructor_of[rft], station_of(s), ph2_dump_node, pr)

    # -------------- open arrival wiring: Source into first step, exits to Sink
    for (eidx, first_step, exits) in open_wiring:
        first_cls = step_class[first_step]
        src_node.setArrival(first_cls, arrival_map[eidx])
        P.set(first_cls, first_cls, src_node, station_of(first_step), 1.0)
        for (s, sig, pr) in exits:
            # Open chains carry no signals, so every exit is an ordinary class.
            ecls = step_class[s]
            P.set(ecls, ecls, station_of(s), snk_node, pr)

    model.link(P)

    # see _kb/04-networkstruct.md (lqn2qn activity-walk mechanics) for rationale
    fcr_list = sorted(fcr_task)
    if fcr_list:
        import numpy as np
        classes = model.getClasses()
        Kc = len(classes)
        Amat = np.zeros((len(fcr_list), Kc))
        bvec = np.zeros(len(fcr_list))
        region_nodes = []
        for ti, tidx in enumerate(fcr_list):
            for i in range(nsteps):
                if tidx not in step_tasks[i]:
                    continue
                cidx = classes.index(step_class[i])
                Amat[ti, cidx] = 1.0
                nd = station_of(i)
                if isinstance(nd, (Queue, Delay)) and nd not in region_nodes:
                    region_nodes.append(nd)
            bvec[ti] = float(lsn.mult[0, tidx])
        if Amat.any() and region_nodes:
            fcr = model.add_region(*region_nodes)
            fcr.set_linear_constraints(Amat, bvec)

    return model


# Alias for snake_case
lqn2qn = LQN2QN
