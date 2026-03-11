"""
Topology validation and parameter extraction for AoI analysis.

This module follows the MATLAB AoI validation and extraction logic used by
the FLD MFQ solver path.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

from ...api.sn import NetworkStruct, NodeType, SchedStrategy


def aoi_is_aoi(sn: NetworkStruct) -> Tuple[bool, Dict[str, Any]]:
    """Validate that `sn` matches the supported AoI topology."""
    info: Dict[str, Any] = {
        'errorMsg': '',
        'sourceIdx': None,
        'queueIdx': None,
        'sinkIdx': None,
        'sourceStation': None,
        'queueStation': None,
        'capacity': None,
        'schedStrategy': None,
        'systemType': '',
        'openClass': None,
    }

    njobs = np.asarray(getattr(sn, 'njobs', np.array([]))).reshape(-1)
    open_classes = np.flatnonzero(np.isinf(njobs))
    if open_classes.size == 0 and getattr(sn, 'nclasses', 0) == 1 and getattr(sn, 'nclosedjobs', 0) == 0:
        open_classes = np.array([0], dtype=int)

    if open_classes.size == 0:
        return _aoi_error(info, 'Not an open model - all classes are closed')
    if open_classes.size > 1:
        return _aoi_error(
            info,
            f"Multiple open classes found ({open_classes.size}) - AoI analysis requires single class",
        )

    nodetype = np.asarray(getattr(sn, 'nodetype', np.array([]))).reshape(-1)
    source_nodes = np.flatnonzero(nodetype == int(NodeType.SOURCE))
    queue_nodes = np.flatnonzero(nodetype == int(NodeType.QUEUE))
    sink_nodes = np.flatnonzero(nodetype == int(NodeType.SINK))

    if source_nodes.size == 0:
        return _aoi_error(info, 'No source node found')
    if source_nodes.size > 1:
        return _aoi_error(info, f'Multiple source nodes found ({source_nodes.size})')
    if sink_nodes.size == 0:
        return _aoi_error(info, 'No sink node found')
    if sink_nodes.size > 1:
        return _aoi_error(info, f'Multiple sink nodes found ({sink_nodes.size})')
    if queue_nodes.size == 0:
        return _aoi_error(info, 'No queue node found')
    if queue_nodes.size > 1:
        return _aoi_error(
            info,
            f"Multiple queue nodes found ({queue_nodes.size}) - AoI analysis supports single queue only",
        )

    node_to_station = np.asarray(getattr(sn, 'nodeToStation', np.array([]))).reshape(-1)
    source_idx = int(source_nodes[0])
    queue_idx = int(queue_nodes[0])
    sink_idx = int(sink_nodes[0])
    if source_idx >= len(node_to_station) or queue_idx >= len(node_to_station):
        return _aoi_error(info, 'Invalid node-to-station mapping for AoI topology')

    source_station = int(node_to_station[source_idx])
    queue_station = int(node_to_station[queue_idx])
    if source_station < 0 or queue_station < 0:
        return _aoi_error(info, 'AoI analysis requires source and queue to be stations')

    info.update({
        'sourceIdx': source_idx,
        'queueIdx': queue_idx,
        'sinkIdx': sink_idx,
        'sourceStation': source_station,
        'queueStation': queue_station,
        'openClass': int(open_classes[0]),
    })

    nservers = np.asarray(getattr(sn, 'nservers', np.array([]))).reshape(-1)
    if nservers.size and queue_station < len(nservers) and int(nservers[queue_station]) != 1:
        return _aoi_error(
            info,
            f'Queue has {int(nservers[queue_station])} servers - AoI analysis requires single server (c=1)',
        )

    cap_value = _station_value(sn, 'cap', queue_station)
    if cap_value is None:
        cap_value = _station_value(sn, 'capacity', queue_station)
    if cap_value is None or np.isinf(cap_value) or cap_value > 2 or cap_value < 1:
        return _aoi_error(
            info,
            f'Queue capacity is {cap_value} - AoI analysis requires capacity 1 (bufferless) or 2 (single-buffer)',
        )

    capacity = int(cap_value)
    sched_strategy = _station_sched(sn, queue_station)
    if sched_strategy not in (SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.LCFSPR):
        return _aoi_error(
            info,
            'Unsupported scheduling strategy - AoI analysis supports FCFS, LCFS, or LCFSPR only',
        )

    if capacity == 2 and not _is_exponential_arrival(sn, source_station, int(open_classes[0])):
        return _aoi_error(
            info,
            'Single-buffer (capacity=2) requires exponential arrivals (Poisson process)',
        )

    rt = np.asarray(getattr(sn, 'rt', np.array([])))
    if rt.size:
        k = int(getattr(sn, 'nclasses', 0))
        rt_idx = queue_station * k + int(open_classes[0])
        if rt_idx < rt.shape[0] and rt_idx < rt.shape[1] and rt[rt_idx, rt_idx] > 0:
            return _aoi_error(info, 'Self-loop detected at queue - violates AoI model assumptions')

    info['capacity'] = capacity
    info['schedStrategy'] = sched_strategy
    info['systemType'] = 'bufferless' if capacity == 1 else 'singlebuffer'
    info['is_valid'] = True
    info['reason'] = ''
    return True, info


def aoi_extract_params(sn: NetworkStruct, aoi_info: Dict[str, Any], options: Any) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Extract the PH/MFQ parameters required by the AoI solvers."""
    from .conversion import aoi_dist2ph

    source_station = int(aoi_info['sourceStation'])
    queue_station = int(aoi_info['queueStation'])
    class_idx = int(aoi_info['openClass'])
    capacity = int(aoi_info['capacity'])
    sched_strategy = aoi_info['schedStrategy']

    service_proc = _get_process(getattr(sn, 'proc', None), queue_station, class_idx)
    if service_proc is None:
        mu = _station_rate(sn, queue_station, class_idx, fallback=1.0)
        service_proc = _exp_dist(mu)
    sigma, S = aoi_dist2ph(service_proc)

    override = getattr(options, 'aoi_preemption', None)
    if capacity == 1:
        arrival_proc = _get_process(getattr(sn, 'proc', None), source_station, class_idx)
        if arrival_proc is None:
            lambda_rate = _station_rate(sn, source_station, class_idx, fallback=0.5)
            arrival_proc = _exp_dist(lambda_rate)
        tau, T = aoi_dist2ph(arrival_proc)

        if override is not None:
            p = float(override)
        elif sched_strategy == SchedStrategy.LCFSPR:
            p = 1.0
        else:
            p = 0.0

        params = {'tau': tau, 'T': T, 'sigma': sigma, 'S': S, 'p': p}
        config = {
            'systemType': 'bufferless',
            'queue_idx': queue_station,
            'source_idx': source_station,
            'preemption_type': 'LCFSPR' if p > 0.5 else 'FCFS',
        }
    else:
        lambda_rate = _station_rate(sn, source_station, class_idx, fallback=0.5)
        if override is not None:
            r = float(override)
        elif sched_strategy in (SchedStrategy.LCFS, SchedStrategy.LCFSPR):
            r = 1.0
        else:
            r = 0.0

        params = {'lambda_rate': lambda_rate, 'sigma': sigma, 'S': S, 'r': r}
        config = {
            'systemType': 'singlebuffer',
            'queue_idx': queue_station,
            'source_idx': source_station,
            'replacement_type': 'replacement' if r > 0.5 else 'FCFS',
        }

    return params, config


def _aoi_error(info: Dict[str, Any], message: str) -> Tuple[bool, Dict[str, Any]]:
    info['errorMsg'] = message
    info['reason'] = message
    info['is_valid'] = False
    return False, info


def _station_value(sn: NetworkStruct, attr: str, station: int) -> Optional[float]:
    value = getattr(sn, attr, None)
    if value is None:
        return None
    arr = np.asarray(value).reshape(-1)
    if station >= len(arr):
        return None
    return float(arr[station])


def _station_rate(sn: NetworkStruct, station: int, class_idx: int, fallback: float) -> float:
    rates = getattr(sn, 'rates', None)
    if rates is not None:
        rates_arr = np.asarray(rates, dtype=float)
        if rates_arr.ndim == 2 and station < rates_arr.shape[0] and class_idx < rates_arr.shape[1]:
            return float(rates_arr[station, class_idx])
    lambda_arr = getattr(sn, 'lambda_arr', None)
    if lambda_arr is not None:
        lam_arr = np.asarray(lambda_arr, dtype=float).reshape(-1)
        if class_idx < len(lam_arr):
            return float(lam_arr[class_idx])
    return float(fallback)


def _station_sched(sn: NetworkStruct, station: int) -> int:
    sched_dict = getattr(sn, 'sched', None)
    if isinstance(sched_dict, dict) and station in sched_dict:
        return int(sched_dict[station])
    schedid = getattr(sn, 'schedid', None)
    if schedid is not None:
        schedid_arr = np.asarray(schedid).reshape(-1)
        if station < len(schedid_arr):
            return int(schedid_arr[station])
    return int(SchedStrategy.FCFS)


def _get_process(proc_container: Any, station: int, class_idx: int) -> Any:
    if proc_container is None:
        return None
    if isinstance(proc_container, dict):
        for key in ((station, class_idx), station):
            if key not in proc_container:
                continue
            value = proc_container[key]
            if isinstance(value, dict):
                if class_idx in value:
                    return value[class_idx]
            elif isinstance(value, (list, tuple)):
                if class_idx < len(value):
                    return value[class_idx]
            else:
                return value
        return None
    if isinstance(proc_container, (list, tuple)) and station < len(proc_container):
        station_entry = proc_container[station]
        if isinstance(station_entry, dict):
            return station_entry.get(class_idx)
        if isinstance(station_entry, (list, tuple)) and class_idx < len(station_entry):
            return station_entry[class_idx]
    return None


def _is_exponential_arrival(sn: NetworkStruct, source_station: int, class_idx: int) -> bool:
    arrival_proc = _get_process(getattr(sn, 'proc', None), source_station, class_idx)
    if arrival_proc is None:
        return True
    if hasattr(arrival_proc, 'D0'):
        d0 = np.asarray(arrival_proc.D0)
        return d0.shape == (1, 1)
    proc_name = type(arrival_proc).__name__
    if proc_name in {'Exp', 'Exponential'}:
        return True
    if hasattr(arrival_proc, 'alpha') and hasattr(arrival_proc, 'A'):
        amat = np.asarray(arrival_proc.A)
        return amat.shape == (1, 1)
    if hasattr(arrival_proc, 'cv'):
        return abs(float(arrival_proc.cv) - 1.0) < 1e-8
    return False


def _exp_dist(rate: float):
    mean = 1.0 / float(rate)
    return type('ExpDist', (), {'mean': mean, 'cv': 1.0})()
