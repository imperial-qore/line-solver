"""Static applicability test for the 'mfq' fluid method and its aliases.

Python twin of the MATLAB ``fluid_mfq_admits`` and ``fluid_is_single_queue``,
kept in one module because the second exists only to serve the first.

The Markovian fluid queue answers ONE open station, Source -> Queue -> Sink.
Both this port and the reference used to fall back to the matrix method from
inside the mfq arm, so a caller's 'mfq' label was answered by a different
algorithm under that label, and the feature set could not say which. The
predicate is stated here instead, so the SAME verdict is available to the
analyzer (which warns and re-enters the matrix method) and to the gate
(``SolverFLD.resolveMethod``, which relabels the pair as 'matrix'). One
predicate, two callers, as in MATLAB.

Three arms, tried in the analyzer's own order:

* the age-of-information shape (``aoi_is_aoi``): one open class, one server, a
  buffer of 1 or 2, FCFS/LCFS/LCFSPR. The finite buffer IS the model there, so
  the third return value lets a caller exempt it from the binding-capacity
  gate, as 'mol' is exempt for the same reason.
* the priority branch, taken when the classes carry distinct priorities: a
  class-independent service rate, a MAP {D0,D1} arrival for every open class,
  at least one of them modulated.
* the plain branch, which reads the arrival and service processes of class 1
  only, so it is stated for ONE open class.

Closed classes are refused by the feature set (``SolverFLD.getMethodFeatureSet``)
rather than here: the plain branch skips them and reports zeros.

NOTE FOR THE PRIORITY ARM. MATLAB and C++ carry a priority solver behind it
(``solver_mfq_prio``, ``fluid_mfq_prio.h``); this port does not, so a model the
priority arm admits is answered by the plain branch. That is what this port did
before the predicate existed too, so the arm is kept faithful to the reference
rather than narrowed here; see ``_kb/07-cross-language-parity.md``.
"""

from typing import Any, Dict, Tuple

import numpy as np

from line_solver.api.sn import NodeType
from line_solver.constants import GlobalConstants


def fluid_is_single_queue(sn) -> Tuple[bool, Dict[str, Any]]:
    """Whether `sn` is the Source -> Queue -> Sink shape 'mfq' is stated for.

    Parameters
    ----------
    sn : NetworkStruct
        Network structure from ``Network.getStruct()``.

    Returns
    -------
    (is_single_queue, info)
        ``info`` carries the node and station indices on success, and
        ``info['errorMsg']`` names the blocking shape otherwise.
    """
    info: Dict[str, Any] = {
        'errorMsg': '',
        'sourceIdx': None,
        'queueIdx': None,
        'sinkIdx': None,
        'sourceStation': None,
        'queueStation': None,
    }

    njobs = np.asarray(getattr(sn, 'njobs', np.array([]))).reshape(-1)
    if njobs.size == 0 or not np.any(np.isinf(njobs)):
        return _shape_error(info, 'Not an open model - all classes are closed')

    nodetype = np.asarray(getattr(sn, 'nodetype', np.array([]))).reshape(-1)
    source_nodes = np.flatnonzero(nodetype == int(NodeType.SOURCE))
    queue_nodes = np.flatnonzero(nodetype == int(NodeType.QUEUE))
    sink_nodes = np.flatnonzero(nodetype == int(NodeType.SINK))

    if source_nodes.size == 0:
        return _shape_error(info, 'No source node found')
    if source_nodes.size > 1:
        return _shape_error(info, 'Multiple source nodes found (%d)' % source_nodes.size)
    if sink_nodes.size == 0:
        return _shape_error(info, 'No sink node found')
    if sink_nodes.size > 1:
        return _shape_error(info, 'Multiple sink nodes found (%d)' % sink_nodes.size)
    if queue_nodes.size == 0:
        return _shape_error(info, 'No queue node found')
    if queue_nodes.size > 1:
        return _shape_error(
            info,
            'Multiple queue nodes found (%d) - MFQ supports single queue only' % queue_nodes.size)

    node_to_station = np.asarray(getattr(sn, 'nodeToStation', np.array([]))).reshape(-1)
    source_idx = int(source_nodes[0])
    queue_idx = int(queue_nodes[0])
    sink_idx = int(sink_nodes[0])
    if source_idx >= node_to_station.size or queue_idx >= node_to_station.size:
        return _shape_error(info, 'Invalid node-to-station mapping for the MFQ topology')

    source_station = int(node_to_station[source_idx])
    queue_station = int(node_to_station[queue_idx])
    if source_station < 0 or queue_station < 0:
        return _shape_error(info, 'MFQ requires source and queue to be stations')

    # c = 1 or c = Inf. The drift of the mfq arm carries one server or infinitely
    # many; min(n,c) is the matrix method's, which is where a c-server model goes.
    nservers = np.asarray(getattr(sn, 'nservers', np.array([]))).reshape(-1)
    if nservers.size and queue_station < nservers.size:
        c = float(nservers[queue_station])
        if c > 1 and np.isfinite(c):
            return _shape_error(
                info,
                'Multi-server queue (c=%g) not supported - MFQ requires c=1 or c=Inf' % c)

    # No self-loop at the queue: the fluid queue is a renewal input to one server
    # and a job fed back to it violates that independence.
    rt = np.asarray(getattr(sn, 'rt', np.array([])))
    if rt.size:
        nclasses = int(getattr(sn, 'nclasses', 0))
        for k in range(nclasses):
            rt_idx = queue_station * nclasses + k
            if rt_idx < rt.shape[0] and rt_idx < rt.shape[1] and rt[rt_idx, rt_idx] > 0:
                return _shape_error(
                    info, 'Self-loop detected at queue - violates independence assumption')

    info.update({
        'sourceIdx': source_idx,
        'queueIdx': queue_idx,
        'sinkIdx': sink_idx,
        'sourceStation': source_station,
        'queueStation': queue_station,
    })
    return True, info


def fluid_mfq_admits(sn) -> Tuple[bool, str, bool]:
    """Whether the 'mfq' method (and its alias 'aoi') may run here.

    Parameters
    ----------
    sn : NetworkStruct
        Network structure from ``Network.getStruct()``.

    Returns
    -------
    (ok, reason, is_aoi)
        ``ok`` is True when 'mfq' may run as itself; ``reason`` names the
        blocking shape otherwise and is empty when ``ok``; ``is_aoi`` is True
        when the age-of-information arm is the one that runs.
    """
    from line_solver.lib.thirdparty.aoi import aoi_is_aoi

    is_aoi, _ = aoi_is_aoi(sn)
    if is_aoi:
        return True, '', True

    ok, info = fluid_is_single_queue(sn)
    if not ok:
        return (False,
                "The 'mfq' method is a single-queue fluid model (Source -> Queue -> Sink, one or "
                "infinitely many servers): %s." % info['errorMsg'],
                False)

    njobs = np.asarray(getattr(sn, 'njobs', np.array([]))).reshape(-1)
    open_classes = np.flatnonzero(np.isinf(njobs))

    classprio = np.asarray(getattr(sn, 'classprio', np.array([]))).reshape(-1)
    # the same test the analyzer switches on
    if classprio.size and np.unique(classprio).size > 1:
        return _admits_priority(sn, info, open_classes)

    if open_classes.size != 1:
        return (False,
                "The 'mfq' method reads the arrival and service processes of one open class and "
                "this model has %d; several classes are served only by its priority branch "
                "(distinct class priorities)." % open_classes.size,
                False)
    return True, '', False


def _admits_priority(sn, info, open_classes) -> Tuple[bool, str, bool]:
    """The priority branch's own three conditions, in the reference's order."""
    qi = info['queueStation']
    si = info['sourceStation']
    rates = np.asarray(getattr(sn, 'rates', np.array([])))
    if open_classes.size == 0 or rates.size == 0 or qi >= rates.shape[0]:
        return (False,
                "The priority branch of 'mfq' needs a finite positive service rate for every "
                "open class.",
                False)
    mu = np.asarray([rates[qi, k] for k in open_classes], dtype=float)
    if not np.all(np.isfinite(mu)) or np.any(mu <= 0):
        return (False,
                "The priority branch of 'mfq' needs a finite positive service rate for every "
                "open class.",
                False)
    if np.any(np.abs(mu - mu[0]) > GlobalConstants.FineTol * max(1.0, mu[0])):
        return (False,
                "The priority branch of 'mfq' is a fluid priority queue drained at one rate, so "
                "it needs a class-independent service rate.",
                False)

    nph = 1
    proc = getattr(sn, 'proc', None)
    for k in open_classes:
        entry = None
        if proc is not None:
            try:
                entry = proc[si][int(k)]
            except (IndexError, KeyError, TypeError):
                entry = None
        if entry is None or len(entry) < 2 or np.asarray(entry[0]).ndim != 2:
            return (False,
                    "The priority branch of 'mfq' needs a MAP {D0,D1} arrival process for "
                    "class %d." % (int(k) + 1),
                    False)
        nph *= int(np.asarray(entry[0]).shape[0])
    if nph < 2:
        return (False,
                "The priority branch of 'mfq' needs at least one Markov-modulated (multi-phase) "
                "arrival process: with exponential arrivals its fluid model degenerates.",
                False)
    return True, '', False


def _shape_error(info: Dict[str, Any], msg: str) -> Tuple[bool, Dict[str, Any]]:
    info['errorMsg'] = msg
    return False, info
