"""
SN utility functions.

Native Python implementations of network structure utility functions
for printing, debugging, and data conversion.

Port from:
    - matlab/src/api/sn/sn_print.m
    - matlab/src/api/sn/sn_print_routing_matrix.m
    - matlab/src/api/sn/sn_refresh_process_fields.m
    - matlab/src/api/sn/sn_rtnodes_to_rtorig.m
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple, List
import sys

from .network_struct import NetworkStruct, NodeType, RoutingStrategy


def sn_print(sn: NetworkStruct, file=None) -> None:
    """
    Print comprehensive information about a NetworkStruct object.

    This function displays all fields, matrices, lists, and maps in a formatted
    manner useful for debugging and inspection of network structures.

    Args:
        sn: Network structure to inspect
        file: Output file (default: sys.stdout)

    References:
        MATLAB: matlab/src/api/sn/sn_print.m
    """
    if file is None:
        file = sys.stdout

    def p(s):
        print(s, file=file)

    def format_matrix(m, name=''):
        """Format a matrix for compact display."""
        if m is None:
            return 'null'
        if isinstance(m, np.ndarray):
            if m.size == 0:
                return '[]'
            if m.ndim == 1:
                items = [str(int(x)) if x == int(x) and not np.isinf(x) else str(x) for x in m]
                return '[' + ' '.join(items) + ']'
            else:
                rows = []
                for i in range(m.shape[0]):
                    items = [str(int(x)) if x == int(x) and not np.isinf(x) else str(x) for x in m[i, :]]
                    rows.append(' '.join(items))
                return '[' + '; '.join(rows) + ']'
        return str(m)

    # Basic integer fields
    p(f'nstations: {sn.nstations}')
    p(f'nstateful: {sn.nstateful}')
    p(f'nnodes: {sn.nnodes}')
    p(f'nclasses: {sn.nclasses}')
    p(f'nclosedjobs: {sn.nclosedjobs}')
    p(f'nchains: {sn.nchains}')

    # Matrix fields
    p(f'refstat: {format_matrix(sn.refstat)}')
    p(f'njobs: {format_matrix(sn.njobs)}')
    p(f'nservers: {format_matrix(sn.nservers)}')
    if sn.connmatrix is not None:
        p(f'connmatrix: {format_matrix(sn.connmatrix)}')
    p(f'scv: {format_matrix(sn.scv)}')

    # Mapping arrays
    p(f'nodeToStateful: {format_matrix(sn.nodeToStateful)}')
    p(f'nodeToStation: {format_matrix(sn.nodeToStation)}')
    p(f'stationToNode: {format_matrix(sn.stationToNode)}')
    p(f'stationToStateful: {format_matrix(sn.stationToStateful)}')
    p(f'statefulToStation: {format_matrix(sn.statefulToStation)}')
    p(f'statefulToNode: {format_matrix(sn.statefulToNode)}')

    # Rate fields
    p(f'rates: {format_matrix(sn.rates)}')
    if sn.classprio is not None:
        p(f'classprio: {format_matrix(sn.classprio)}')
    if sn.phases is not None:
        p(f'phases: {format_matrix(sn.phases)}')

    # Node type list
    if sn.nodetype is not None:
        nodetype_names = [NodeType.toText(nt) for nt in sn.nodetype]
        p(f'nodetype: [{", ".join(nodetype_names)}]')

    # Class names
    if sn.classnames is not None:
        if isinstance(sn.classnames, list):
            classnames_str = ', '.join(f'"{n}"' for n in sn.classnames)
            p(f'classnames: [{classnames_str}]')
        else:
            p(f'classnames: ["{sn.classnames}"]')

    # Node names
    if sn.nodenames is not None:
        if isinstance(sn.nodenames, list):
            nodenames_str = ', '.join(f'"{n}"' for n in sn.nodenames)
            p(f'nodenames: [{nodenames_str}]')

    # Routing tables
    if sn.rt is not None:
        p(f'rt: {format_matrix(sn.rt)}')
    if sn.rtnodes is not None:
        p(f'rtnodes: {format_matrix(sn.rtnodes)}')

    # Visit ratios
    if sn.visits:
        p('visits: {')
        for chain_id, visits in sn.visits.items():
            p(f'  {chain_id}: {format_matrix(visits)}')
        p('}')

    # Node visits
    if sn.nodevisits:
        p('nodevisits: {')
        for chain_id, nodevisits in sn.nodevisits.items():
            p(f'  {chain_id}: {format_matrix(nodevisits)}')
        p('}')

    # Chain info
    if sn.inchain:
        p('inchain: {')
        for chain_id, classes in sn.inchain.items():
            if isinstance(classes, np.ndarray):
                classes_str = ', '.join(str(int(c)) for c in classes.flatten())
            elif isinstance(classes, list):
                classes_str = ', '.join(str(c) for c in classes)
            else:
                classes_str = str(classes)
            p(f'  {chain_id}: [{classes_str}]')
        p('}')


def sn_print_routing_matrix(
    sn: NetworkStruct,
    onlyclass: Optional[Any] = None,
    file=None
) -> None:
    """
    Print the routing matrix of the network.

    This function displays the routing probabilities between nodes and classes
    in a human-readable format.

    Args:
        sn: Network structure
        onlyclass: Optional filter for a specific class (object with 'name' attribute)
        file: Output file (default: sys.stdout)

    References:
        MATLAB: matlab/src/api/sn/sn_print_routing_matrix.m
    """
    if file is None:
        file = sys.stdout

    node_names = sn.nodenames if sn.nodenames else [f'Node{i}' for i in range(sn.nnodes)]
    classnames = sn.classnames if sn.classnames else [f'Class{i}' for i in range(sn.nclasses)]
    rtnodes = sn.rtnodes
    nnodes = sn.nnodes
    nclasses = sn.nclasses

    if rtnodes is None:
        print("No routing matrix available.", file=file)
        return

    for i in range(nnodes):
        for r in range(nclasses):
            for j in range(nnodes):
                for s in range(nclasses):
                    rt_idx_src = i * nclasses + r
                    rt_idx_dst = j * nclasses + s

                    if rt_idx_src < rtnodes.shape[0] and rt_idx_dst < rtnodes.shape[1]:
                        prob = rtnodes[rt_idx_src, rt_idx_dst]
                        if prob > 0:
                            # Skip sinks
                            if sn.nodetype is not None and i < len(sn.nodetype) and sn.nodetype[i] == NodeType.SINK:
                                continue

                            # Check for Cache (state-dependent)
                            if sn.nodetype is not None and i < len(sn.nodetype) and sn.nodetype[i] == NodeType.CACHE:
                                pr_str = 'state-dependent'
                            elif sn.routing is not None and i < sn.routing.shape[0] and r < sn.routing.shape[1]:
                                if sn.routing[i, r] == RoutingStrategy.DISABLED:
                                    continue
                                pr_str = f'{prob:.6f}'
                            else:
                                pr_str = f'{prob:.6f}'

                            # Apply class filter if specified
                            if onlyclass is None:
                                print(f'{node_names[i]} [{classnames[r]}] => {node_names[j]} [{classnames[s]}] : Pr={pr_str}', file=file)
                            else:
                                class_name = getattr(onlyclass, 'name', str(onlyclass))
                                if classnames[r].lower() == class_name.lower() or classnames[s].lower() == class_name.lower():
                                    print(f'{node_names[i]} [{classnames[r]}] => {node_names[j]} [{classnames[s]}] : Pr={pr_str}', file=file)


def sn_refresh_process_fields(
    sn: NetworkStruct,
    station_idx: int,
    class_idx: int
) -> NetworkStruct:
    """
    Refresh process fields based on rate and SCV values.

    Updates mu, phi, proc, pie, phases based on current rate and SCV values.
    - SCV = 1.0: Exponential (1 phase)
    - SCV < 1.0: Erlang approximation
    - SCV > 1.0: Hyperexponential(2) approximation

    Args:
        sn: Network structure (modified in place)
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)

    Returns:
        Modified network structure

    References:
        MATLAB: matlab/src/api/sn/sn_refresh_process_fields.m
    """
    if sn.rates is None or station_idx >= sn.rates.shape[0] or class_idx >= sn.rates.shape[1]:
        return sn

    rate = sn.rates[station_idx, class_idx]
    scv = sn.scv[station_idx, class_idx] if sn.scv is not None else 1.0

    # Skip if rate is invalid
    if np.isnan(rate) or rate <= 0 or np.isinf(rate):
        return sn

    mean = 1.0 / rate

    # Determine number of phases and create MAP representation
    if np.isnan(scv) or abs(scv - 1.0) < 1e-10:
        # Exponential
        n_phases = 1
        D0 = np.array([[-rate]])
        D1 = np.array([[rate]])
    elif scv < 1.0:
        # Erlang: k = ceil(1/scv)
        k = max(1, int(np.ceil(1.0 / scv)))
        n_phases = k
        # Erlang-k phase rate
        phase_rate = k * rate
        D0 = np.zeros((k, k))
        for i in range(k - 1):
            D0[i, i] = -phase_rate
            D0[i, i + 1] = phase_rate
        D0[k - 1, k - 1] = -phase_rate
        D1 = np.zeros((k, k))
        D1[k - 1, 0] = phase_rate
    else:
        # Hyperexponential (scv > 1)
        n_phases = 2
        # Fit H2 to match mean and scv
        cv2 = scv
        if cv2 <= 1:
            cv2 = 1.01

        # Two-moment matching for H2
        # Using balanced means approach
        p = 0.5 * (1 + np.sqrt((cv2 - 1) / (cv2 + 1)))
        p = min(max(p, 0.01), 0.99)

        mu1 = 2 * p / mean
        mu2 = 2 * (1 - p) / mean

        D0 = np.array([[-mu1, 0], [0, -mu2]])
        D1 = np.array([[mu1, 0], [0, mu2]])

    # Update phases
    if sn.phases is None:
        sn.phases = np.ones((sn.nstations, sn.nclasses))
    if station_idx < sn.phases.shape[0] and class_idx < sn.phases.shape[1]:
        sn.phases[station_idx, class_idx] = n_phases

    # Update phasessz
    if sn.phasessz is None:
        sn.phasessz = np.ones((sn.nstations, sn.nclasses))
    if station_idx < sn.phasessz.shape[0] and class_idx < sn.phasessz.shape[1]:
        sn.phasessz[station_idx, class_idx] = max(n_phases, 1)

    # Recompute phaseshift for this station
    if sn.phaseshift is None:
        sn.phaseshift = np.zeros((sn.nstations, sn.nclasses + 1))
    cum_sum = 0
    sn.phaseshift[station_idx, 0] = 0
    for c in range(sn.nclasses):
        if sn.phasessz is not None:
            cum_sum += sn.phasessz[station_idx, c]
        if c + 1 < sn.phaseshift.shape[1]:
            sn.phaseshift[station_idx, c + 1] = cum_sum

    # Update mu (rates from -diag(D0))
    if sn.mu is None:
        sn.mu = {}
    mu_vec = -np.diag(D0)
    sn.mu[(station_idx, class_idx)] = mu_vec

    # Update phi (completion probabilities)
    if sn.phi is None:
        sn.phi = {}
    phi_vec = np.zeros(n_phases)
    for i in range(n_phases):
        d1_row_sum = np.sum(D1[i, :])
        d0_diag = -D0[i, i]
        if d0_diag != 0:
            phi_vec[i] = d1_row_sum / d0_diag
    sn.phi[(station_idx, class_idx)] = phi_vec

    # Update pie (initial phase distribution)
    if sn.pie is None:
        sn.pie = {}
    pie_vec = np.zeros(n_phases)
    pie_vec[0] = 1.0  # Start in first phase
    sn.pie[(station_idx, class_idx)] = pie_vec

    # Update proc (MAP representation)
    if sn.proc is None:
        sn.proc = {}
    sn.proc[(station_idx, class_idx)] = [D0, D1]

    return sn


def sn_rtnodes_to_rtorig(sn: NetworkStruct) -> Tuple[Dict, np.ndarray]:
    """
    Convert node routing matrix to the original routing matrix format.

    This function converts the node-level routing matrix to the original
    routing matrix format, excluding class-switching nodes.

    Args:
        sn: Network structure

    Returns:
        Tuple of (rtorigcell, rtorig) where:
            rtorigcell: Dictionary representation {(r,s): ndarray}
            rtorig: Sparse/dense matrix representation

    References:
        MATLAB: matlab/src/api/sn/sn_rtnodes_to_rtorig.m
    """
    K = sn.nclasses
    rtnodes = sn.rtnodes

    if rtnodes is None:
        return {}, np.array([])

    # Find where class-switching nodes start
    csshift = sn.nnodes
    if sn.nodenames is not None:
        for ind in range(sn.nnodes):
            if ind < len(sn.nodenames) and sn.nodenames[ind].startswith('CS_'):
                csshift = ind
                break

    # Build column indices to keep (exclude CS nodes)
    col_to_keep = []
    for ind in range(csshift):
        for k in range(K):
            col_to_keep.append(ind * K + k)

    if len(col_to_keep) == 0:
        return {}, np.array([])

    # Perform stochastic complementation
    rtorig = _dtmc_stochcomp(rtnodes, col_to_keep)

    # Replace NaNs with 0
    rtorig = np.nan_to_num(rtorig, nan=0.0)

    # Build cell representation
    rtorigcell: Dict[Tuple[int, int], np.ndarray] = {}
    for r in range(K):
        for s in range(K):
            rtorigcell[(r, s)] = np.zeros((csshift, csshift))

    for ind in range(csshift):
        if sn.nodetype is not None and ind < len(sn.nodetype) and sn.nodetype[ind] != NodeType.SINK:
            for jnd in range(csshift):
                for r in range(K):
                    for s in range(K):
                        src_idx = ind * K + r
                        dst_idx = jnd * K + s
                        if src_idx < rtorig.shape[0] and dst_idx < rtorig.shape[1]:
                            rtorigcell[(r, s)][ind, jnd] = rtorig[src_idx, dst_idx]

    return rtorigcell, rtorig


def _dtmc_stochcomp(P: np.ndarray, keep_states: List[int]) -> np.ndarray:
    """
    Perform stochastic complementation on a transition matrix.

    Removes transient states by computing the stochastic complement.

    Args:
        P: Transition probability matrix
        keep_states: Indices of states to keep

    Returns:
        Reduced transition probability matrix
    """
    n = P.shape[0]
    keep_set = set(keep_states)
    remove_states = [i for i in range(n) if i not in keep_set]

    if len(remove_states) == 0:
        # Nothing to remove
        return P[np.ix_(keep_states, keep_states)]

    if len(keep_states) == 0:
        return np.array([])

    # Partition the matrix
    # P = [[Q_AA, Q_AB], [Q_BA, Q_BB]]
    # where A = keep, B = remove
    Q_AA = P[np.ix_(keep_states, keep_states)]
    Q_AB = P[np.ix_(keep_states, remove_states)]
    Q_BA = P[np.ix_(remove_states, keep_states)]
    Q_BB = P[np.ix_(remove_states, remove_states)]

    # Stochastic complement: P_A = Q_AA + Q_AB * (I - Q_BB)^(-1) * Q_BA
    I_BB = np.eye(len(remove_states))
    try:
        inv_term = np.linalg.inv(I_BB - Q_BB)
        P_A = Q_AA + Q_AB @ inv_term @ Q_BA
    except np.linalg.LinAlgError:
        # If inversion fails, just use the direct submatrix
        P_A = Q_AA

    return P_A
