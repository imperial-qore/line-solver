"""
Utility functions to convert between JPype wrapper structures and native Python structures.

This module provides conversion functions to allow native Python solvers to work
with models created through the JPype wrapper.
"""

import numpy as np
from typing import Dict, Optional

from ..api.sn import NetworkStruct as NativeNetworkStruct
from ..api.sn import SchedStrategy as NativeSchedStrategy
from ..api.sn import NodeType as NativeNodeType


# Mapping from string sched names to native SchedStrategy enum values
SCHED_NAME_TO_NATIVE = {
    'FCFS': NativeSchedStrategy.FCFS,
    'LCFS': NativeSchedStrategy.LCFS,
    'LCFSPR': NativeSchedStrategy.LCFSPR,
    'PS': NativeSchedStrategy.PS,
    'DPS': NativeSchedStrategy.DPS,
    'GPS': NativeSchedStrategy.GPS,
    'INF': NativeSchedStrategy.INF,
    'RAND': NativeSchedStrategy.RAND,
    'HOL': NativeSchedStrategy.HOL,
    'SEPT': NativeSchedStrategy.SEPT,
    'LEPT': NativeSchedStrategy.LEPT,
    'SIRO': NativeSchedStrategy.SIRO,
    'SJF': NativeSchedStrategy.SJF,
    'LJF': NativeSchedStrategy.LJF,
    'POLLING': NativeSchedStrategy.POLLING,
    'EXT': NativeSchedStrategy.EXT,
}

# Mapping from string node type names to native NodeType enum values
NODE_TYPE_NAME_TO_NATIVE = {
    'SOURCE': NativeNodeType.SOURCE,
    'SINK': NativeNodeType.SINK,
    'QUEUE': NativeNodeType.QUEUE,
    'DELAY': NativeNodeType.DELAY,
    'ROUTER': NativeNodeType.ROUTER,
    'FORK': NativeNodeType.FORK,
    'JOIN': NativeNodeType.JOIN,
    'CACHE': NativeNodeType.CACHE,
    'LOGGER': NativeNodeType.LOGGER,
    'CLASSSWITCH': NativeNodeType.CLASSSWITCH,
    'PLACE': NativeNodeType.PLACE,
    'TRANSITION': NativeNodeType.TRANSITION,
    'FINITE_CAPACITY_REGION': NativeNodeType.FINITE_CAPACITY_REGION,
}


def get_native_sched_strategy(sched_name: str) -> int:
    """
    Convert scheduling strategy string name to native enum value.

    Args:
        sched_name: Scheduling strategy name (e.g., 'FCFS', 'PS', 'INF')

    Returns:
        Native SchedStrategy enum value
    """
    return SCHED_NAME_TO_NATIVE.get(sched_name, NativeSchedStrategy.FCFS)


def get_native_node_type(node_type) -> int:
    """
    Convert node type to native enum value.

    Args:
        node_type: Node type (can be enum, string, or int)

    Returns:
        Native NodeType enum value
    """
    if isinstance(node_type, int):
        return node_type
    elif hasattr(node_type, 'name'):
        # Case-insensitive lookup
        name = node_type.name.upper()
        return NODE_TYPE_NAME_TO_NATIVE.get(name, NativeNodeType.QUEUE)
    elif isinstance(node_type, str):
        return NODE_TYPE_NAME_TO_NATIVE.get(node_type.upper(), NativeNodeType.QUEUE)
    return NativeNodeType.QUEUE


def wrapper_sn_to_native(sn) -> NativeNetworkStruct:
    """
    Convert a JPype wrapper NetworkStruct to native Python NetworkStruct.

    This function takes a NetworkStruct from the JPype wrapper (line_solver.lang.NetworkStruct)
    and converts it to a native Python NetworkStruct (line_solver.api.sn.NetworkStruct)
    that can be used by native Python solvers.

    Args:
        sn: NetworkStruct from JPype wrapper (line_solver.lang.NetworkStruct)
            Can also be a model object with getStruct() method

    Returns:
        NativeNetworkStruct suitable for native solvers
    """
    # If passed a model, get its struct first
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct(force=True)

    # Build sched dict from ndarray
    sched_dict = {}
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist in range(sn.nstations):
            if ist < len(sn.sched):
                sched_name = sn.sched[ist]
                if isinstance(sched_name, str):
                    sched_dict[ist] = get_native_sched_strategy(sched_name)
                elif sched_name is None:
                    sched_dict[ist] = NativeSchedStrategy.FCFS
                elif hasattr(sched_name, 'value'):
                    # It's an enum, use its value
                    sched_dict[ist] = int(sched_name.value)
                elif hasattr(sched_name, 'name'):
                    # It's an enum, convert by name
                    sched_dict[ist] = get_native_sched_strategy(sched_name.name)
                else:
                    sched_dict[ist] = int(sched_name)

    # Convert node types
    nodetype_list = []
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for nt in sn.nodetype:
            nodetype_list.append(get_native_node_type(nt))

    # Convert visits dict
    visits_dict = {}
    if hasattr(sn, 'visits') and sn.visits is not None:
        for chain_id, visit_matrix in sn.visits.items():
            if visit_matrix is not None:
                visits_dict[chain_id] = np.asarray(visit_matrix, dtype=np.float64)

    # Convert inchain dict
    inchain_dict = {}
    if hasattr(sn, 'inchain') and sn.inchain is not None:
        for chain_id, class_indices in sn.inchain.items():
            if class_indices is not None:
                inchain_dict[chain_id] = np.asarray(class_indices, dtype=np.int64)

    # Helper to safely convert arrays
    def safe_array(arr, dtype=np.float64):
        if arr is None:
            return np.array([])
        return np.asarray(arr, dtype=dtype)

    # Helper to safely convert names
    def safe_names(names_list):
        if names_list is None:
            return []
        result = []
        for name in names_list:
            if name is None:
                result.append('')
            elif isinstance(name, str):
                result.append(name)
            elif hasattr(name, 'getName'):
                # Object with getName() method
                result.append(str(name.getName()))
            elif hasattr(name, 'name'):
                # Python object with name attribute
                result.append(str(name.name))
            else:
                # Fallback: convert to string
                result.append(str(name))
        return result

    # Create native NetworkStruct
    native_sn = NativeNetworkStruct(
        nstations=int(sn.nstations),
        nstateful=int(sn.nstateful) if hasattr(sn, 'nstateful') else int(sn.nstations),
        nnodes=int(sn.nnodes) if hasattr(sn, 'nnodes') else int(sn.nstations),
        nclasses=int(sn.nclasses),
        nchains=int(sn.nchains) if hasattr(sn, 'nchains') else 1,
        nclosedjobs=int(sn.nclosedjobs) if hasattr(sn, 'nclosedjobs') else 0,

        # Population and capacity
        njobs=safe_array(sn.njobs),
        nservers=safe_array(sn.nservers),
        cap=safe_array(sn.cap) if hasattr(sn, 'cap') and sn.cap is not None else None,
        classcap=safe_array(sn.classcap) if hasattr(sn, 'classcap') and sn.classcap is not None else None,

        # Service parameters
        rates=safe_array(sn.rates),
        scv=safe_array(sn.scv),
        phases=safe_array(sn.phases) if hasattr(sn, 'phases') and sn.phases is not None else None,

        # Chains
        visits=visits_dict,
        inchain=inchain_dict,
        chains=safe_array(sn.chains) if hasattr(sn, 'chains') else np.array([]),
        refstat=safe_array(sn.refstat, dtype=np.int64).flatten(),
        refclass=safe_array(sn.refclass) if hasattr(sn, 'refclass') else np.array([]),

        # Scheduling and routing
        sched=sched_dict,
        schedparam=safe_array(sn.schedparam) if hasattr(sn, 'schedparam') and sn.schedparam is not None else None,
        rt=safe_array(sn.rt) if hasattr(sn, 'rt') and sn.rt is not None else None,
        rtnodes=safe_array(sn.rtnodes) if hasattr(sn, 'rtnodes') and sn.rtnodes is not None else None,

        # Node classification
        nodetype=nodetype_list,
        isstation=safe_array(sn.isstation) if hasattr(sn, 'isstation') else np.array([]),
        isstateful=safe_array(sn.isstateful) if hasattr(sn, 'isstateful') else np.array([]),

        # Mappings (flatten to 1D arrays)
        nodeToStation=safe_array(sn.nodeToStation, dtype=np.int64).flatten() if hasattr(sn, 'nodeToStation') else np.array([]),
        nodeToStateful=safe_array(sn.nodeToStateful, dtype=np.int64).flatten() if hasattr(sn, 'nodeToStateful') else np.array([]),
        stationToNode=safe_array(sn.stationToNode, dtype=np.int64).flatten() if hasattr(sn, 'stationToNode') else np.array([]),
        stationToStateful=safe_array(sn.stationToStateful, dtype=np.int64).flatten() if hasattr(sn, 'stationToStateful') else np.array([]),
        statefulToNode=safe_array(sn.statefulToNode, dtype=np.int64).flatten() if hasattr(sn, 'statefulToNode') else np.array([]),
        statefulToStation=safe_array(sn.statefulToStation, dtype=np.int64).flatten() if hasattr(sn, 'statefulToStation') else np.array([]),

        # Load-dependent scaling
        lldscaling=safe_array(sn.lldscaling) if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else None,

        # Class properties
        classprio=safe_array(sn.classprio) if hasattr(sn, 'classprio') and sn.classprio is not None else None,

        # Connectivity
        connmatrix=safe_array(sn.connmatrix) if hasattr(sn, 'connmatrix') else None,

        # Convert names to strings if needed
        nodenames=safe_names(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else [],
        classnames=safe_names(sn.classnames) if hasattr(sn, 'classnames') and sn.classnames else [],
    )

    # Copy process parameters (proc and procid) if available
    if hasattr(sn, 'proc') and sn.proc is not None:
        native_sn.proc = _convert_proc(sn.proc, native_sn.nstations, native_sn.nclasses)

    if hasattr(sn, 'procid') and sn.procid is not None:
        native_sn.procid = _convert_procid(sn.procid, native_sn.nstations, native_sn.nclasses)

    return native_sn


def _convert_proc(proc, nstations: int, nclasses: int):
    """
    Convert proc field from wrapper format to native format.

    The proc field contains distribution parameters (D0, D1 matrices for MAP/MMPP2,
    alpha/T for PH, etc.) indexed by [station][class].

    Args:
        proc: Process parameters from wrapper NetworkStruct
        nstations: Number of stations
        nclasses: Number of classes

    Returns:
        Native proc structure (list of lists)
    """
    if proc is None:
        return None

    # Initialize native proc structure
    native_proc = [[None for _ in range(nclasses)] for _ in range(nstations)]

    try:
        # Handle cell array from Java/MATLAB (indexed as proc[i][j] or proc{i}{j})
        for i in range(nstations):
            for j in range(nclasses):
                try:
                    # Try various access patterns
                    proc_ij = None

                    # Direct indexing (list of lists or numpy array)
                    if hasattr(proc, '__getitem__'):
                        if i < len(proc) and proc[i] is not None:
                            if hasattr(proc[i], '__getitem__') and j < len(proc[i]):
                                proc_ij = proc[i][j]

                    if proc_ij is not None:
                        # Convert to native format
                        if isinstance(proc_ij, (list, tuple)) and len(proc_ij) >= 2:
                            # [D0, D1] format for MAP/MMPP2
                            D0 = np.asarray(proc_ij[0], dtype=np.float64) if proc_ij[0] is not None else None
                            D1 = np.asarray(proc_ij[1], dtype=np.float64) if proc_ij[1] is not None else None
                            if D0 is not None and D1 is not None:
                                native_proc[i][j] = [D0, D1]
                        elif hasattr(proc_ij, 'items'):
                            # Dict format
                            native_proc[i][j] = dict(proc_ij)
                        else:
                            # Other format - try to convert
                            native_proc[i][j] = proc_ij

                except (IndexError, TypeError, KeyError):
                    continue

    except Exception:
        # If conversion fails, return None
        return None

    return native_proc


def _convert_procid(procid, nstations: int, nclasses: int):
    """
    Convert procid field from wrapper format to native format.

    The procid field contains ProcessType values indexed by [station, class].

    Args:
        procid: Process type IDs from wrapper NetworkStruct
        nstations: Number of stations
        nclasses: Number of classes

    Returns:
        Native procid array (numpy array of ProcessType)
    """
    from ..constants import ProcessType

    if procid is None:
        return None

    # Initialize with default EXP
    native_procid = np.empty((nstations, nclasses), dtype=object)
    native_procid.fill(ProcessType.EXP)

    try:
        procid_arr = np.asarray(procid)

        for i in range(min(nstations, procid_arr.shape[0])):
            for j in range(min(nclasses, procid_arr.shape[1] if len(procid_arr.shape) > 1 else 1)):
                try:
                    val = procid_arr[i, j] if len(procid_arr.shape) > 1 else procid_arr[i]

                    # Convert to ProcessType
                    if isinstance(val, ProcessType):
                        native_procid[i, j] = val
                    elif hasattr(val, 'name'):
                        # Enum from Java
                        proc_type = ProcessType.fromString(val.name)
                        if proc_type is not None:
                            native_procid[i, j] = proc_type
                    elif isinstance(val, str):
                        proc_type = ProcessType.fromString(val)
                        if proc_type is not None:
                            native_procid[i, j] = proc_type
                    elif isinstance(val, (int, np.integer)):
                        # Try to map integer to ProcessType
                        for pt in ProcessType:
                            if pt.value == int(val):
                                native_procid[i, j] = pt
                                break

                except (IndexError, TypeError, KeyError):
                    continue

    except Exception:
        # If conversion fails, return default
        pass

    return native_procid


def native_model_to_struct(model) -> NativeNetworkStruct:
    """
    Convert pure Python model to NetworkStruct

    This function takes various types of Python model representations and
    converts them to a NativeNetworkStruct that can be used by native solvers.

    Supports:
    - lang.Network objects (pure Python)
    - Dict-based model specifications
    - Pre-built NetworkStruct objects

    Args:
        model: Model object. Can be:
            - lang.Network with get_struct() or refresh_struct()
            - Dict-based model specification
            - NativeNetworkStruct (returned as-is)
            - JPype wrapper Network (falls back to wrapper_sn_to_native)

    Returns:
        NativeNetworkStruct suitable for native solvers

    Raises:
        TypeError: If model cannot be converted to NetworkStruct
    """
    # Case 1: Already a NetworkStruct
    if isinstance(model, NativeNetworkStruct):
        return model

    # Case 2: Native Network with _sn attribute
    if hasattr(model, '_sn') and model._sn is not None:
        if isinstance(model._sn, NativeNetworkStruct):
            return model._sn
        # Might be a different struct type, try to convert
        return wrapper_sn_to_native(model._sn)

    # Case 3: Native Network with refresh_struct()
    if hasattr(model, 'refresh_struct'):
        model.refresh_struct()
        if hasattr(model, '_sn') and model._sn is not None:
            if isinstance(model._sn, NativeNetworkStruct):
                return model._sn
            return wrapper_sn_to_native(model._sn)

    # Case 4: Native Network with get_struct()
    if hasattr(model, 'get_struct'):
        sn = model.get_struct()
        if isinstance(sn, NativeNetworkStruct):
            return sn
        return wrapper_sn_to_native(sn)

    # Case 5: JPype wrapper with getStruct() method
    if hasattr(model, 'getStruct'):
        sn = model.getStruct(force=True)
        return wrapper_sn_to_native(sn)

    # Case 6: JPype wrapper with obj attribute
    if hasattr(model, 'obj'):
        try:
            sn = model.getStruct(force=True)
            return wrapper_sn_to_native(sn)
        except Exception:
            pass

    # Case 7: Dict-based model
    if isinstance(model, dict):
        return _build_struct_from_dict(model)

    # Case 8: Already a struct-like object with required attributes
    if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
        return wrapper_sn_to_native(model)

    raise TypeError(f"Cannot convert {type(model)} to NetworkStruct. "
                    "Expected lang.Network, dict, or NetworkStruct-like object.")


def _build_struct_from_dict(model_dict: dict) -> NativeNetworkStruct:
    """
    Build NetworkStruct from dict specification.

    Dict keys should match NetworkStruct field names:
    - nstations, nclasses (required)
    - rates, njobs, nservers, visits, etc. (optional)

    Args:
        model_dict: Dictionary with network parameters

    Returns:
        NativeNetworkStruct

    Raises:
        ValueError: If required keys are missing
    """
    # Check required keys
    if 'nstations' not in model_dict:
        raise ValueError("Dict model must include 'nstations'")
    if 'nclasses' not in model_dict:
        raise ValueError("Dict model must include 'nclasses'")

    nstations = int(model_dict['nstations'])
    nclasses = int(model_dict['nclasses'])

    # Helper to safely convert arrays
    def safe_array(key, default_shape=None, dtype=np.float64):
        if key in model_dict and model_dict[key] is not None:
            return np.asarray(model_dict[key], dtype=dtype)
        if default_shape is not None:
            return np.zeros(default_shape, dtype=dtype)
        return np.array([], dtype=dtype)

    # Build sched dict from dict or array
    sched_dict = {}
    if 'sched' in model_dict and model_dict['sched'] is not None:
        sched_val = model_dict['sched']
        if isinstance(sched_val, dict):
            for k, v in sched_val.items():
                if isinstance(v, str):
                    sched_dict[int(k)] = get_native_sched_strategy(v)
                else:
                    sched_dict[int(k)] = int(v) if v is not None else NativeSchedStrategy.FCFS
        elif hasattr(sched_val, '__iter__'):
            for i, v in enumerate(sched_val):
                if isinstance(v, str):
                    sched_dict[i] = get_native_sched_strategy(v)
                else:
                    sched_dict[i] = int(v) if v is not None else NativeSchedStrategy.FCFS

    # Build visits dict
    visits_dict = {}
    if 'visits' in model_dict and model_dict['visits'] is not None:
        visits_val = model_dict['visits']
        if isinstance(visits_val, dict):
            for k, v in visits_val.items():
                if v is not None:
                    visits_dict[int(k)] = np.asarray(v, dtype=np.float64)

    # Build inchain dict
    inchain_dict = {}
    if 'inchain' in model_dict and model_dict['inchain'] is not None:
        inchain_val = model_dict['inchain']
        if isinstance(inchain_val, dict):
            for k, v in inchain_val.items():
                if v is not None:
                    inchain_dict[int(k)] = np.asarray(v, dtype=np.int64)

    # Convert node types if present
    nodetype_list = []
    if 'nodetype' in model_dict and model_dict['nodetype'] is not None:
        for nt in model_dict['nodetype']:
            nodetype_list.append(get_native_node_type(nt))

    # Create native NetworkStruct
    native_sn = NativeNetworkStruct(
        nstations=nstations,
        nstateful=int(model_dict.get('nstateful', nstations)),
        nnodes=int(model_dict.get('nnodes', nstations)),
        nclasses=nclasses,
        nchains=int(model_dict.get('nchains', 1)),
        nclosedjobs=int(model_dict.get('nclosedjobs', 0)),

        # Population and capacity
        njobs=safe_array('njobs', (nclasses, 1)),
        nservers=safe_array('nservers', (nstations, 1)),
        cap=safe_array('cap') if 'cap' in model_dict else None,
        classcap=safe_array('classcap') if 'classcap' in model_dict else None,

        # Service parameters
        rates=safe_array('rates', (nstations, nclasses)),
        scv=safe_array('scv', (nstations, nclasses)),
        phases=safe_array('phases') if 'phases' in model_dict else None,

        # Chains
        visits=visits_dict,
        inchain=inchain_dict,
        chains=safe_array('chains'),
        refstat=safe_array('refstat', dtype=np.int64).flatten(),
        refclass=safe_array('refclass'),

        # Scheduling and routing
        sched=sched_dict,
        schedparam=safe_array('schedparam') if 'schedparam' in model_dict else None,
        rt=safe_array('rt') if 'rt' in model_dict else None,
        rtnodes=safe_array('rtnodes') if 'rtnodes' in model_dict else None,

        # Node classification
        nodetype=nodetype_list,
        isstation=safe_array('isstation'),
        isstateful=safe_array('isstateful'),

        # Mappings
        nodeToStation=safe_array('nodeToStation', dtype=np.int64).flatten(),
        nodeToStateful=safe_array('nodeToStateful', dtype=np.int64).flatten(),
        stationToNode=safe_array('stationToNode', dtype=np.int64).flatten(),
        stationToStateful=safe_array('stationToStateful', dtype=np.int64).flatten(),
        statefulToNode=safe_array('statefulToNode', dtype=np.int64).flatten(),
        statefulToStation=safe_array('statefulToStation', dtype=np.int64).flatten(),

        # Load-dependent scaling
        lldscaling=safe_array('lldscaling') if 'lldscaling' in model_dict else None,

        # Class properties
        classprio=safe_array('classprio') if 'classprio' in model_dict else None,

        # Connectivity
        connmatrix=safe_array('connmatrix') if 'connmatrix' in model_dict else None,

        # Names
        nodenames=list(model_dict.get('nodenames', [])),
        classnames=list(model_dict.get('classnames', [])),
    )

    return native_sn


__all__ = [
    'wrapper_sn_to_native',
    'native_model_to_struct',
    'get_native_sched_strategy',
    'get_native_node_type',
]
