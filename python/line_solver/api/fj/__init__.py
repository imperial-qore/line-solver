"""
Fork-Join (FJ) Topology Analysis.

Native Python implementations for analyzing Fork-Join queueing systems,
including topology validation, parameter extraction, distribution mapping,
and synchronization map construction.

Key functions:
    fj_isfj: Check if network has valid Fork-Join topology for FJ_codes analysis
    fj_extract_params: Extract Fork-Join parameters from network
    fj_dist2fj: Convert LINE MAP distribution to FJ_codes format
    sn_build_fj_sync_map: Build fork-join synchronization map for MAM solver

References:
    Original MATLAB: matlab/src/api/fj/fj_*.m
    Qiu, Perez, Harrison, "Beyond the Mean in Fork-Join Queues", 2015
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any, List
from dataclasses import dataclass, field
import warnings


@dataclass
class FjParams:
    """Fork-Join network parameters."""
    K: int  # Number of parallel queues
    lambda_val: float  # Arrival rate
    mu: np.ndarray  # Service rates at each queue
    cv: np.ndarray  # Coefficient of variation at each queue
    dist_type: str  # Distribution type ('Exp', 'Erlang', 'HyperExp')


@dataclass
class FjInfo:
    """Fork-Join topology information returned by fj_isfj.

    Matches MATLAB struct returned by fj_isfj.m.
    """
    forkIdx: Optional[int] = None   # Fork node index (0-based)
    joinIdx: Optional[int] = None   # Join node index (0-based)
    queueIdx: Optional[List[int]] = field(default_factory=list)  # Queue node indices between Fork and Join
    K: int = 0                      # Number of parallel queues
    errorMsg: str = ''              # Error message if not valid FJ
    arrival: Optional[list] = None  # Arrival process parameters (one per class)
    service: Optional[list] = None  # Service process parameters (one per class)
    nClasses: int = 0               # Number of classes


@dataclass
class FjSyncMap:
    """Fork-Join synchronization map for MAM solver.

    Matches MATLAB struct returned by sn_build_fj_sync_map.m.
    """
    nodeSync: np.ndarray = None     # nodeSync[joinIdx, srcIdx] = groupId
    forkOfGroup: list = field(default_factory=list)  # forkOfGroup[groupId] = forkIdx
    joinOfGroup: list = field(default_factory=list)  # joinOfGroup[groupId] = joinIdx
    nGroups: int = 0                # Total number of sync groups


def fj_isfj(sn) -> Tuple[bool, FjInfo]:
    """
    Check if network has valid Fork-Join topology for FJ_codes analysis.

    Validates that the network structure matches the requirements for FJ_codes:
    - Single Fork-Join pair
    - K parallel queues between Fork and Join
    - Homogeneous service distributions across parallel queues
    - Supported distributions (Exp, HyperExp(2), Erlang(2), MAP(2))
    - Open classes only
    - FCFS or PS scheduling

    Args:
        sn: Network structure (NetworkStruct)

    Returns:
        Tuple of (isFJ, fjInfo) where:
            isFJ: True if network is valid FJ topology for FJ_codes
            fjInfo: FjInfo with fields forkIdx, joinIdx, queueIdx, K, errorMsg

    References:
        Original MATLAB: matlab/src/api/fj/fj_isfj.m
        Z. Qiu, J.F. Perez, and P. Harrison, "Beyond the Mean in Fork-Join
        Queues: Efficient Approximation for Response-Time Tails", IFIP
        Performance 2015.
    """
    from ..sn.network_struct import NodeType, SchedStrategy
    from ..sn.predicates import sn_is_open_model, sn_has_fork_join

    fjInfo = FjInfo()

    # Check if model has open classes only
    if not sn_is_open_model(sn):
        fjInfo.errorMsg = 'FJ_codes only supports open queueing models.'
        return False, fjInfo

    # Check if network has fork-join
    if not sn_has_fork_join(sn):
        fjInfo.errorMsg = 'Network does not contain Fork-Join structure.'
        return False, fjInfo

    # Find Fork and Join nodes
    nodetype = sn.nodetype
    if nodetype is None:
        fjInfo.errorMsg = 'Network has no nodetype information.'
        return False, fjInfo

    fork_indices = []
    join_indices = []
    for i in range(len(nodetype)):
        if nodetype[i] == NodeType.FORK:
            fork_indices.append(i)
        elif nodetype[i] == NodeType.JOIN:
            join_indices.append(i)

    if len(fork_indices) == 0 or len(join_indices) == 0:
        fjInfo.errorMsg = 'Network must contain both Fork and Join nodes.'
        return False, fjInfo

    # FJ_codes supports single Fork-Join pair
    if len(fork_indices) > 1:
        fjInfo.errorMsg = 'FJ_codes only supports a single Fork-Join pair. Found multiple Fork nodes.'
        return False, fjInfo

    if len(join_indices) > 1:
        fjInfo.errorMsg = 'FJ_codes only supports a single Fork-Join pair. Found multiple Join nodes.'
        return False, fjInfo

    forkIdx = fork_indices[0]
    joinIdx = join_indices[0]

    # Check if Fork and Join are paired using sn.fj matrix
    if sn.fj is None or sn.fj[forkIdx, joinIdx] == 0:
        fjInfo.errorMsg = (
            f'Fork node {forkIdx} and Join node {joinIdx} are not paired.'
        )
        return False, fjInfo

    fjInfo.forkIdx = forkIdx
    fjInfo.joinIdx = joinIdx

    # Find queues between Fork and Join
    # These are nodes that receive routing from Fork and route to Join
    queueIdx = []
    K = sn.nclasses
    rtnodes = sn.rtnodes
    for i in range(sn.nnodes):
        if nodetype[i] == NodeType.QUEUE:
            hasForkInput = False
            hasJoinOutput = False

            # Check if Fork routes to this queue (any class)
            for k in range(K):
                fork_row = forkIdx * K + k
                node_col = i * K + k
                if rtnodes is not None and fork_row < rtnodes.shape[0] and node_col < rtnodes.shape[1]:
                    if rtnodes[fork_row, node_col] > 0:
                        hasForkInput = True

                node_row = i * K + k
                join_col = joinIdx * K + k
                if rtnodes is not None and node_row < rtnodes.shape[0] and join_col < rtnodes.shape[1]:
                    if rtnodes[node_row, join_col] > 0:
                        hasJoinOutput = True

            if hasForkInput and hasJoinOutput:
                queueIdx.append(i)

    if len(queueIdx) == 0:
        fjInfo.errorMsg = 'No Queue nodes found between Fork and Join.'
        return False, fjInfo

    num_queues = len(queueIdx)
    fjInfo.queueIdx = queueIdx
    fjInfo.K = num_queues

    # Validate homogeneous service distributions across parallel queues
    FineTol = 1e-8
    for r in range(sn.nclasses):
        firstQueueIdx = queueIdx[0]
        firstQueueStat = sn.nodeToStation[firstQueueIdx]
        if firstQueueStat < 0:
            fjInfo.errorMsg = f'Queue {firstQueueIdx} is not a station.'
            return False, fjInfo

        firstPH = None
        if sn.proc is not None:
            firstPH = sn.proc[firstQueueStat][r]

        if firstPH is None or (isinstance(firstPH, list) and len(firstPH) > 0 and
                                isinstance(firstPH[0], np.ndarray) and np.isnan(firstPH[0].flat[0])):
            fjInfo.errorMsg = (
                f'Queue {firstQueueIdx} has no valid service distribution for class {r}.'
            )
            return False, fjInfo

        # Check all other queues have the same distribution
        for ki in range(1, num_queues):
            queueIdx_k = queueIdx[ki]
            queueStat_k = sn.nodeToStation[queueIdx_k]
            if queueStat_k < 0:
                fjInfo.errorMsg = f'Queue {queueIdx_k} is not a station.'
                return False, fjInfo

            ph_k = None
            if sn.proc is not None:
                ph_k = sn.proc[queueStat_k][r]

            if ph_k is None or (isinstance(ph_k, list) and len(ph_k) > 0 and
                                 isinstance(ph_k[0], np.ndarray) and np.isnan(ph_k[0].flat[0])):
                fjInfo.errorMsg = (
                    f'Queue {queueIdx_k} has no valid service distribution for class {r}.'
                )
                return False, fjInfo

            # Compare PH representations (must be identical)
            if isinstance(firstPH, list) and isinstance(ph_k, list):
                if len(firstPH) != len(ph_k):
                    fjInfo.errorMsg = (
                        f'Queues have heterogeneous service distributions for class {r}. '
                        f'FJ_codes requires homogeneous servers.'
                    )
                    return False, fjInfo

                for mat_idx in range(min(len(firstPH), len(ph_k))):
                    if isinstance(firstPH[mat_idx], np.ndarray) and isinstance(ph_k[mat_idx], np.ndarray):
                        if firstPH[mat_idx].shape != ph_k[mat_idx].shape:
                            fjInfo.errorMsg = (
                                f'Queues have heterogeneous service distributions for class {r}. '
                                f'FJ_codes requires homogeneous servers.'
                            )
                            return False, fjInfo
                        if np.max(np.abs(firstPH[mat_idx] - ph_k[mat_idx])) > FineTol:
                            fjInfo.errorMsg = (
                                f'Queues have heterogeneous service distributions for class {r}. '
                                f'FJ_codes requires homogeneous servers.'
                            )
                            return False, fjInfo

    # Validate supported scheduling strategies (FCFS or PS)
    for ki in range(num_queues):
        queueStat = sn.nodeToStation[queueIdx[ki]]
        if sn.sched is not None:
            sched_val = sn.sched[queueStat]
            if sched_val != SchedStrategy.FCFS and sched_val != SchedStrategy.PS:
                fjInfo.errorMsg = (
                    f'Queue {queueIdx[ki]} has unsupported scheduling strategy. '
                    f'FJ_codes supports FCFS or PS only.'
                )
                return False, fjInfo

    # All validations passed
    return True, fjInfo


def fj_extract_params_detailed(sn, fjInfo):
    """
    Extract FJ_codes parameters from LINE network structure.

    Extracts arrival process, service process, and K value from a validated
    Fork-Join network structure for use with FJ_codes analysis.

    Args:
        sn: Network structure (after conversion)
        fjInfo: FjInfo from fj_isfj (with forkIdx, joinIdx, queueIdx, K)

    Returns:
        Tuple of (arrival, service, K, fjInfo) where:
            arrival: list of dicts (one per class) with fields:
                lambda_val, lambda0, lambda1, ma, Ia, ArrChoice
            service: list of dicts (one per class) with fields:
                mu, ST, St, tau_st, SerChoice
            K: Number of parallel queues
            fjInfo: Updated FjInfo

    References:
        Original MATLAB: matlab/src/api/fj/fj_extract_params.m
    """
    from ..sn.network_struct import NodeType

    forkIdx = fjInfo.forkIdx
    joinIdx = fjInfo.joinIdx
    queueIdx = fjInfo.queueIdx
    K = fjInfo.K

    # Find Source node to extract arrival process
    nodetype = sn.nodetype
    sourceIdx = None
    for i in range(len(nodetype)):
        if nodetype[i] == NodeType.SOURCE:
            sourceIdx = i
            break

    if sourceIdx is None:
        raise ValueError('No Source node found in network.')

    sourceStat = sn.nodeToStation[sourceIdx]

    # Initialize lists for multi-class support
    nClasses = sn.nclasses
    arrival = [None] * nClasses
    service = [None] * nClasses

    # Extract parameters for each class
    for r in range(nClasses):
        # === ARRIVAL PROCESS ===
        arrivalMAP = None
        if sn.proc is not None:
            arrivalMAP = sn.proc[sourceStat][r]

        if arrivalMAP is None:
            raise ValueError(f'Source has no valid arrival process for class {r}.')

        arrival[r] = fj_dist2fj_map(arrivalMAP, 'arrival', sn, sourceStat, r)

        # === SERVICE PROCESS ===
        firstQueue = queueIdx[0]
        firstQueueStat = sn.nodeToStation[firstQueue]
        serviceMAP = None
        if sn.proc is not None:
            serviceMAP = sn.proc[firstQueueStat][r]

        if serviceMAP is None:
            raise ValueError(
                f'Queue {firstQueue} has no valid service process for class {r}.'
            )

        service[r] = fj_dist2fj_map(serviceMAP, 'service', sn, firstQueueStat, r)

        # === STABILITY CHECK ===
        arr_rate = arrival[r].get('lambda_val', 0)
        svc_rate = service[r].get('mu', float('inf'))
        if arr_rate >= svc_rate:
            warnings.warn(
                f'Class {r} may be unstable: arrival rate ({arr_rate:.4f}) >= '
                f'service rate ({svc_rate:.4f}). FJ_codes requires stable '
                f'systems (rho = lambda/mu < 1).'
            )

    # Store distribution info in fjInfo
    fjInfo.arrival = arrival
    fjInfo.service = service
    fjInfo.nClasses = nClasses

    return arrival, service, K, fjInfo


def fj_extract_params(sn) -> Optional[FjParams]:
    """
    Extract Fork-Join parameters from network structure.

    Simplified interface that extracts basic parameters for FJ analysis.
    For full MAP-based extraction, use fj_extract_params_detailed.

    Args:
        sn: Network structure

    Returns:
        FjParams if valid FJ network, None otherwise

    References:
        Original MATLAB: matlab/src/api/fj/fj_extract_params.m
    """
    from ..sn.network_struct import NodeType

    isFJ, fjInfo = fj_isfj(sn)
    if not isFJ:
        return None

    # Find source and extract arrival rate
    nodetype = sn.nodetype
    source_idx = None
    for i in range(len(nodetype)):
        if nodetype[i] == NodeType.SOURCE:
            source_idx = i
            break

    if source_idx is None:
        return None

    # Get arrival rate
    lambda_val = 0.0
    if hasattr(sn, 'rates') and sn.rates is not None:
        sourceStat = sn.nodeToStation[source_idx]
        if sourceStat >= 0:
            lambda_val = np.sum(sn.rates[sourceStat, :])

    # Find queues and extract service parameters
    K = fjInfo.K
    queueIdx = fjInfo.queueIdx

    if K == 0:
        return None

    mu = np.zeros(K)
    cv = np.ones(K)  # Default CV = 1 (exponential)
    dist_type = 'Exp'

    for i, q_idx in enumerate(queueIdx):
        ist = sn.nodeToStation[q_idx]
        if ist >= 0 and hasattr(sn, 'rates') and sn.rates is not None:
            for r in range(sn.nclasses):
                rate = sn.rates[ist, r]
                if not np.isnan(rate) and rate > 0:
                    mu[i] = rate
                    # Get SCV for CV
                    if hasattr(sn, 'scv') and sn.scv is not None:
                        scv_val = sn.scv[ist, r]
                        if not np.isnan(scv_val) and scv_val >= 0:
                            cv[i] = np.sqrt(scv_val)
                    break

    return FjParams(
        K=K,
        lambda_val=lambda_val,
        mu=mu,
        cv=cv,
        dist_type=dist_type
    )


def fj_dist2fj_map(lineMAP, dist_type, sn, ist, r):
    """
    Convert LINE PH/MAP distribution to FJ_codes format.

    Converts LINE's MAP representation {D0, D1} (stored as a list [D0, D1])
    to FJ_codes format.

    For arrivals: returns dict with lambda_val, lambda0, lambda1, ma, Ia, ArrChoice
    For service: returns dict with mu, ST, St, tau_st, SerChoice

    Supported distributions:
    - Exponential (Exp)
    - 2-phase Hyperexponential (HyperExp with 2 phases)
    - 2-phase Erlang (Erlang with 2 phases)
    - 2-phase MAP (MAP with 2 states)

    Args:
        lineMAP: LINE MAP representation [D0, D1] (list of 2 numpy arrays)
        dist_type: 'arrival' or 'service'
        sn: Network structure (for getting rate and process type)
        ist: Station index
        r: Class index

    Returns:
        Dictionary with FJ_codes compatible parameters

    References:
        Original MATLAB: matlab/src/api/fj/fj_dist2fj.m
        Z. Qiu, J.F. Perez, and P. Harrison, "Beyond the Mean in Fork-Join
        Queues", IFIP Performance 2015.
    """
    result = {}

    # Get process type
    procType = None
    if hasattr(sn, 'procid') and sn.procid is not None:
        procType = sn.procid[ist, r]

    # Extract D0 and D1 from LINE MAP format
    D0 = np.atleast_2d(np.asarray(lineMAP[0], dtype=float))
    D1 = np.atleast_2d(np.asarray(lineMAP[1], dtype=float))
    nPhases = D0.shape[0]

    if nPhases > 2:
        raise ValueError(
            f'FJ_codes only supports distributions with at most 2 phases. '
            f'Found {nPhases} phases.'
        )

    # Compute mean arrival/service rate using map_lambda
    try:
        from ..mam.map_analysis import map_lambda, map_pie
        lambda_rate = map_lambda(D0, D1)
    except (ImportError, Exception):
        # Fallback: compute from generator
        D = D0 + D1
        n = D.shape[0]
        # Stationary distribution of D
        if n == 1:
            lambda_rate = float(D1[0, 0])
        else:
            # Solve pi * D = 0, sum(pi) = 1
            A = D.T.copy()
            A[-1, :] = 1.0
            b = np.zeros(n)
            b[-1] = 1.0
            try:
                pi = np.linalg.solve(A, b)
            except np.linalg.LinAlgError:
                pi = np.ones(n) / n
            lambda_rate = float(pi @ D1 @ np.ones(n))

    if dist_type == 'arrival':
        # === ARRIVAL PROCESS ===
        result['lambda_val'] = lambda_rate
        result['lambda0'] = D0.copy()
        result['lambda1'] = D1.copy()
        result['ma'] = nPhases
        result['Ia'] = np.eye(nPhases)

        # Determine arrival choice for FJ_codes
        if procType is not None:
            try:
                from ...constants import ProcessType
                if procType == ProcessType.EXP or procType == ProcessType.EXP.value:
                    result['ArrChoice'] = 1  # Exp
                elif (procType == ProcessType.HYPEREXP or procType == ProcessType.HYPEREXP.value) and nPhases == 2:
                    result['ArrChoice'] = 2  # HE2
                elif (procType == ProcessType.ERLANG or procType == ProcessType.ERLANG.value) and nPhases == 2:
                    result['ArrChoice'] = 3  # ER2
                elif (procType == ProcessType.MAP or procType == ProcessType.MAP.value) and nPhases == 2:
                    result['ArrChoice'] = 4  # MAP2
                else:
                    result['ArrChoice'] = 0  # Unknown
            except ImportError:
                result['ArrChoice'] = 0
        else:
            # Infer from nPhases
            if nPhases == 1:
                result['ArrChoice'] = 1  # Exp
            elif nPhases == 2:
                result['ArrChoice'] = 4  # MAP2 (most general 2-phase)
            else:
                result['ArrChoice'] = 0

    else:
        # === SERVICE PROCESS ===
        # Compute initial probability vector (stationary distribution of MAP)
        try:
            from ..mam.map_analysis import map_pie
            pie = map_pie(D0, D1)
        except (ImportError, Exception):
            n = D0.shape[0]
            if n == 1:
                pie = np.array([1.0])
            else:
                D = D0 + D1
                A = D.T.copy()
                A[-1, :] = 1.0
                b = np.zeros(n)
                b[-1] = 1.0
                try:
                    pie = np.linalg.solve(A, b)
                except np.linalg.LinAlgError:
                    pie = np.ones(n) / n

        result['mu'] = lambda_rate  # Mean service rate
        result['ST'] = D0.copy()    # Rate matrix
        result['St'] = -D0 @ np.ones((nPhases, 1))  # Exit vector
        result['tau_st'] = pie.flatten()  # Initial probability (row vector)

        # Determine service choice
        if procType is not None:
            try:
                from ...constants import ProcessType
                if procType == ProcessType.EXP or procType == ProcessType.EXP.value:
                    result['SerChoice'] = 1  # Exp
                elif (procType == ProcessType.HYPEREXP or procType == ProcessType.HYPEREXP.value) and nPhases == 2:
                    result['SerChoice'] = 2  # HE2
                elif (procType == ProcessType.ERLANG or procType == ProcessType.ERLANG.value) and nPhases == 2:
                    result['SerChoice'] = 3  # ER2
                else:
                    result['SerChoice'] = 0  # Unknown
            except ImportError:
                result['SerChoice'] = 0
        else:
            if nPhases == 1:
                result['SerChoice'] = 1  # Exp
            elif nPhases == 2:
                result['SerChoice'] = 2  # H2 (default 2-phase)
            else:
                result['SerChoice'] = 0

    return result


def fj_dist2fj(dist_type_or_map, params=None):
    """
    Convert a distribution to FJ-compatible format.

    This function supports two calling conventions:

    1. Simple string-based: fj_dist2fj('exp', {'rate': 2.0})
    2. MAP-based: fj_dist2fj_map(lineMAP, dist_type, sn, ist, r)
       (use fj_dist2fj_map directly for this)

    Supported distributions (string-based):
    - Exp: Exponential (rate)
    - Erlang: Erlang-k (k, rate)
    - HyperExp: 2-phase hyperexponential (p, mu1, mu2)
    - MAP: 2-state MAP (D0, D1)
    - Det: Deterministic (value)

    Args:
        dist_type_or_map: Distribution type string or MAP representation
        params: Distribution parameters (dict)

    Returns:
        Dictionary with FJ-compatible parameters:
            - type: Distribution type
            - mean: Mean value
            - cv: Coefficient of variation
            - phases: Number of phases (if applicable)
            - representation: Matrix representation (if applicable)

    References:
        Original MATLAB: matlab/src/api/fj/fj_dist2fj.m
    """
    dist_type = dist_type_or_map
    if params is None:
        params = {}

    result = {
        'type': dist_type,
        'mean': 0.0,
        'cv': 1.0,
        'phases': 1,
        'representation': None
    }

    if dist_type.lower() in ['exp', 'exponential']:
        rate = params.get('rate', params.get('mu', 1.0))
        result['mean'] = 1.0 / rate
        result['cv'] = 1.0
        result['phases'] = 1

    elif dist_type.lower() in ['erlang', 'erlangk']:
        k = params.get('k', params.get('phases', 2))
        rate = params.get('rate', params.get('mu', 1.0))
        result['mean'] = k / rate
        result['cv'] = 1.0 / np.sqrt(k)
        result['phases'] = k

    elif dist_type.lower() in ['hyperexp', 'hyperexponential', 'h2']:
        p = params.get('p', 0.5)
        mu1 = params.get('mu1', params.get('rate1', 1.0))
        mu2 = params.get('mu2', params.get('rate2', 1.0))

        mean1 = 1.0 / mu1
        mean2 = 1.0 / mu2
        result['mean'] = p * mean1 + (1 - p) * mean2

        # Second moment
        m2 = 2 * (p / mu1 ** 2 + (1 - p) / mu2 ** 2)
        var = m2 - result['mean'] ** 2
        result['cv'] = np.sqrt(var) / result['mean'] if result['mean'] > 0 else 1.0
        result['phases'] = 2

    elif dist_type.lower() in ['map', 'map2']:
        D0 = np.asarray(params.get('D0', [[-1]]))
        D1 = np.asarray(params.get('D1', [[1]]))

        # Compute MAP mean using map_lambda if available
        try:
            from ..mam.map_analysis import map_lambda
            lambda_val = map_lambda(D0, D1)
        except (ImportError, Exception):
            D = D0 + D1
            n = D.shape[0]
            pi = np.ones(n) / n
            for _ in range(100):
                pi_new = pi @ np.linalg.matrix_power(np.eye(n) + D / 100, 100)
                pi_new /= np.sum(pi_new)
                if np.linalg.norm(pi_new - pi) < 1e-12:
                    break
                pi = pi_new
            lambda_val = pi @ D1 @ np.ones(n)

        result['mean'] = 1.0 / lambda_val if lambda_val > 0 else 0
        result['representation'] = {'D0': D0, 'D1': D1}
        result['phases'] = D0.shape[0]

    elif dist_type.lower() in ['det', 'deterministic', 'constant']:
        value = params.get('value', params.get('mean', 1.0))
        result['mean'] = value
        result['cv'] = 0.0
        result['phases'] = float('inf')  # Infinite phases for deterministic

    else:
        # Unknown distribution - use provided mean and CV if available
        result['mean'] = params.get('mean', 1.0)
        result['cv'] = params.get('cv', 1.0)

    return result


def sn_build_fj_sync_map(sn) -> FjSyncMap:
    """
    Build a fork-join synchronization map from LINE's sn structure.

    For each (Fork, Join) pair, identifies which source nodes feed into the
    join and must be synchronized using mmap_max.

    Args:
        sn: Network structure (NetworkStruct)

    Returns:
        FjSyncMap with fields:
            nodeSync: ndarray [nnodes x nnodes] where
                nodeSync[joinIdx, srcIdx] = groupId
                groupId > 0 means srcIdx belongs to sync group at joinIdx
                groupId == 0 means srcIdx is an independent (non-synced) flow
            forkOfGroup: list where forkOfGroup[groupId] = forkIdx
            joinOfGroup: list where joinOfGroup[groupId] = joinIdx
            nGroups: total number of sync groups

    References:
        Original MATLAB: matlab/src/api/fj/sn_build_fj_sync_map.m
    """
    I = sn.nnodes
    K = sn.nclasses

    # nodeSync is indexed by (destination node, source node)
    nodeSync = np.zeros((I, I), dtype=int)

    groupId = 0
    forkOfGroup = []
    joinOfGroup = []

    # Iterate over all (Fork, Join) pairs from sn.fj
    if sn.fj is None:
        return FjSyncMap(
            nodeSync=nodeSync,
            forkOfGroup=forkOfGroup,
            joinOfGroup=joinOfGroup,
            nGroups=0
        )

    # Find fork indices that have at least one join paired
    fork_indices = []
    for i in range(sn.fj.shape[0]):
        if np.any(sn.fj[i, :] > 0):
            fork_indices.append(i)

    rtnodes = sn.rtnodes

    for forkIdx in fork_indices:
        join_indices = np.where(sn.fj[forkIdx, :] > 0)[0]
        if len(join_indices) == 0:
            continue

        for jnd in join_indices:
            groupId += 1
            forkOfGroup.append(forkIdx)
            joinOfGroup.append(jnd)

            # Find nodes on the parallel path between this fork and join.
            # A node is on the parallel path if:
            #   1. The fork routes to it (for at least one class), AND
            #   2. It routes to the join (for at least one class)
            for ind in range(I):
                if ind == forkIdx or ind == jnd:
                    continue

                forkRoutesToNode = False
                nodeRoutesToJoin = False

                for k in range(K):
                    fork_row = forkIdx * K + k
                    node_col = ind * K + k
                    if (rtnodes is not None and
                            fork_row < rtnodes.shape[0] and
                            node_col < rtnodes.shape[1]):
                        if rtnodes[fork_row, node_col] > 0:
                            forkRoutesToNode = True

                    node_row = ind * K + k
                    join_col = jnd * K + k
                    if (rtnodes is not None and
                            node_row < rtnodes.shape[0] and
                            join_col < rtnodes.shape[1]):
                        if rtnodes[node_row, join_col] > 0:
                            nodeRoutesToJoin = True

                if forkRoutesToNode and nodeRoutesToJoin:
                    nodeSync[jnd, ind] = groupId

    return FjSyncMap(
        nodeSync=nodeSync,
        forkOfGroup=forkOfGroup,
        joinOfGroup=joinOfGroup,
        nGroups=groupId
    )


from .analytical import (
    # Basic functions
    fj_harmonic,
    fj_synch_delay,
    # Response time approximations
    fj_respt_2way,
    fj_respt_nt,
    fj_respt_vm,
    fj_respt_varki,
    fj_rmax,
    fj_bounds,
    fj_rmax_evd,
    fj_rmax_erlang,
    # Expected maximum functions
    fj_xmax_exp,
    fj_xmax_2,
    fj_xmax_erlang,
    fj_xmax_hyperexp,
    fj_xmax_approx,
    fj_xmax_emma,
    fj_xmax_normal,
    fj_xmax_pareto,
    # Order statistics and bounds
    fj_order_stat,
    fj_gk_bound,
    fj_char_max,
    fj_quantile,
    fj_sm_tput,
)

__all__ = [
    # Topology functions
    'FjParams',
    'FjInfo',
    'FjSyncMap',
    'fj_isfj',
    'fj_extract_params',
    'fj_extract_params_detailed',
    'fj_dist2fj',
    'fj_dist2fj_map',
    'sn_build_fj_sync_map',
    # Basic functions
    'fj_harmonic',
    'fj_synch_delay',
    # Response time approximations
    'fj_respt_2way',
    'fj_respt_nt',
    'fj_respt_vm',
    'fj_respt_varki',
    'fj_rmax',
    'fj_bounds',
    'fj_rmax_evd',
    'fj_rmax_erlang',
    # Expected maximum functions
    'fj_xmax_exp',
    'fj_xmax_2',
    'fj_xmax_erlang',
    'fj_xmax_hyperexp',
    'fj_xmax_approx',
    'fj_xmax_emma',
    'fj_xmax_normal',
    'fj_xmax_pareto',
    # Order statistics and bounds
    'fj_order_stat',
    'fj_gk_bound',
    'fj_char_max',
    'fj_quantile',
    'fj_sm_tput',
]
