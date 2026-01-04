"""
QNS Solver handler - Native Python implementation.

Calls the external qnsolver tool via subprocess, bypassing Java/JPype.

Port from:

"""

import numpy as np
import subprocess
import tempfile
import os
import shutil
import platform
from dataclasses import dataclass
from typing import Optional, Tuple, List
from xml.etree import ElementTree as ET
from xml.dom import minidom
import time

from ...sn import (
    NetworkStruct,
    NodeType,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
)


@dataclass
class SolverQNSOptions:
    """Options for QNS solver."""
    method: str = 'default'
    multiserver: str = 'default'
    verbose: bool = False


@dataclass
class SolverQNSReturn:
    """
    Result of QNS solver handler.

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
        it: Number of iterations
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
    method: str = "default"
    it: int = 0


def is_qns_available() -> bool:
    """Check if the qnsolver command is available in PATH."""
    try:
        result = subprocess.run(
            ['qnsolver', '-h'],
            capture_output=True,
            timeout=5
        )
        # qnsolver exists if it returns without error or gives help output
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return False


def _has_multi_server(sn: NetworkStruct) -> bool:
    """Check if the model has multi-server stations."""
    if sn.nservers is None:
        return False
    for i in range(sn.nstations):
        servers = sn.nservers[i] if len(sn.nservers.shape) == 1 else sn.nservers[i, 0]
        if servers > 1 and not np.isinf(servers):
            return True
    return False


def _write_jmva_file(sn: NetworkStruct, model_path: str, options: SolverQNSOptions) -> None:
    """
    Write the network model to JMVA XML format.

    Port from: SolverJMT.writeJMVA() in Java
    """
    # Get chain-level demands
    demands = sn_get_demands_chain(sn)

    M = sn.nstations
    C = sn.nchains

    # Create root element
    root = ET.Element('model')
    root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    root.set('xsi:noNamespaceSchemaLocation', 'JMTmodel.xsd')

    # Algorithm type element
    alg_type = ET.SubElement(root, 'algType')
    alg_type.set('name', 'MVA')
    alg_type.set('tolerance', '1.0E-7')
    alg_type.set('maxSamples', '10000')

    # Parameters element
    params = ET.SubElement(root, 'parameters')

    # Classes element
    classes = ET.SubElement(params, 'classes')
    classes.set('number', str(C))

    # Check if chains contain open or closed classes
    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(sn.nclasses)

    for c in range(C):
        inchain = sn.inchain.get(c, np.array([]))
        if len(inchain) == 0:
            continue

        # Sum njobs for all classes in this chain
        sum_njobs = sum(njobs[int(k)] for k in inchain if int(k) < len(njobs))

        if np.isinf(sum_njobs):
            # Open class - calculate arrival rate
            rate_sum = 0.0
            if sn.rates is not None and sn.nodetype is not None:
                for i in range(len(sn.nodetype)):
                    if sn.nodetype[i] == NodeType.Source:
                        for k in inchain:
                            k_int = int(k)
                            if k_int < sn.rates.shape[1]:
                                val = sn.rates[i, k_int] if len(sn.rates.shape) > 1 else sn.rates[k_int]
                                if not np.isnan(val):
                                    rate_sum += val
            class_elem = ET.SubElement(classes, 'openclass')
            class_elem.set('rate', str(rate_sum))
            class_elem.set('name', f'Chain{c+1:02d}')
        else:
            # Closed class
            class_elem = ET.SubElement(classes, 'closedclass')
            class_elem.set('population', str(int(demands.Nchain[c])))
            class_elem.set('name', f'Chain{c+1:02d}')

    # Stations element
    # Count non-source stations
    source_count = sum(1 for nt in sn.nodetype if nt == NodeType.Source)
    num_stations = M - source_count if sn.nodetype is not None else M

    stations = ET.SubElement(params, 'stations')
    stations.set('number', str(num_stations))

    # Create station elements
    for i in range(M):
        node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
        node_type = sn.nodetype[node_idx] if sn.nodetype is not None else NodeType.Queue

        if node_type == NodeType.Source:
            continue

        station_name = sn.nodenames[node_idx] if sn.nodenames else f'Station{i+1}'

        if node_type == NodeType.Delay:
            stat_elem = ET.SubElement(stations, 'delaystation')
            stat_elem.set('name', station_name)
        elif node_type == NodeType.Queue:
            nservers = sn.nservers[i] if sn.nservers is not None else 1
            if hasattr(nservers, '__len__'):
                nservers = nservers[0] if len(nservers) > 0 else 1

            if nservers == 1 or np.isinf(nservers):
                stat_elem = ET.SubElement(stations, 'listation')
            else:
                stat_elem = ET.SubElement(stations, 'ldstation')
            stat_elem.set('name', station_name)
            stat_elem.set('servers', '1')
        else:
            continue

        # Service times
        srv_times = ET.SubElement(stat_elem, 'servicetimes')
        for c in range(C):
            srv_time = ET.SubElement(srv_times, 'servicetime')
            srv_time.set('customerclass', f'Chain{c+1:02d}')
            st_val = demands.STchain[i, c] if demands.STchain is not None else 0.0
            srv_time.text = str(st_val)

        # Visits
        visits = ET.SubElement(stat_elem, 'visits')
        for c in range(C):
            visit = ET.SubElement(visits, 'visit')
            visit.set('customerclass', f'Chain{c+1:02d}')
            if demands.STchain[i, c] > 0 and demands.Lchain[i, c] > 0:
                visit_val = demands.Lchain[i, c] / demands.STchain[i, c]
            else:
                visit_val = 0.0
            visit.text = str(visit_val)

    # Reference stations
    ref_stations = ET.SubElement(params, 'ReferenceStation')
    ref_stations.set('number', str(C))

    refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(C, dtype=int)
    for c in range(C):
        inchain = sn.inchain.get(c, np.array([]))
        if len(inchain) > 0:
            ref_class_idx = int(inchain[0])
            if ref_class_idx < len(refstat):
                ref_stat_idx = int(refstat[ref_class_idx])
            else:
                ref_stat_idx = 0
        else:
            ref_stat_idx = 0

        node_idx = int(sn.stationToNode[ref_stat_idx]) if sn.stationToNode is not None else ref_stat_idx
        ref_name = sn.nodenames[node_idx] if sn.nodenames and node_idx < len(sn.nodenames) else f'Station{ref_stat_idx+1}'

        class_ref = ET.SubElement(ref_stations, 'Class')
        class_ref.set('name', f'Chain{c+1}')
        class_ref.set('refStation', ref_name)

    # Algorithm parameters
    alg_params = ET.SubElement(root, 'algParams')
    alg_params.append(alg_type)
    compare_algs = ET.SubElement(alg_params, 'compareAlgs')
    compare_algs.set('value', 'false')

    # Write to file with pretty formatting
    xml_str = ET.tostring(root, encoding='unicode')
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')

    # Remove extra blank lines
    lines = [line for line in pretty_xml.split('\n') if line.strip()]
    pretty_xml = '\n'.join(lines)

    with open(model_path, 'w') as f:
        f.write(pretty_xml)


def _parse_results(result_path: str, sn: NetworkStruct, nchains: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse qnsolver output file.

    Returns:
        Tuple of (Uchain, Qchain, Wchain, Tchain)
    """
    M = sn.nstations

    Uchain = np.zeros((M, nchains))
    Qchain = np.zeros((M, nchains))
    Wchain = np.zeros((M, nchains))
    Tchain = np.zeros((M, nchains))

    if not os.path.exists(result_path):
        raise RuntimeError(f"QNS result file not found: {result_path}")

    with open(result_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or '$' in line or line.startswith('#'):
                continue

            parts = line.replace(' ', '').split(',')
            if len(parts) < 5:
                continue

            stat_name = parts[0]

            # Find station index
            station_idx = -1
            for i in range(M):
                node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
                node_name = sn.nodenames[node_idx] if sn.nodenames and node_idx < len(sn.nodenames) else f'Station{i+1}'
                if node_name == stat_name:
                    station_idx = i
                    break

            if station_idx < 0:
                continue

            # Parse values based on output format
            if nchains == 1:
                # Single class format: Station, Q, R, U, X
                try:
                    Qchain[station_idx, 0] = float(parts[1]) if len(parts) > 1 else 0.0
                    Wchain[station_idx, 0] = float(parts[2]) if len(parts) > 2 else 0.0
                    Uchain[station_idx, 0] = float(parts[3]) if len(parts) > 3 else 0.0
                    Tchain[station_idx, 0] = float(parts[4]) if len(parts) > 4 else 0.0
                except (ValueError, IndexError):
                    pass
            else:
                # Multi-class format with aggregate values
                # Format: Station, $Q(Chain01), $Q(Chain02), $Q, $R(Chain01), ...
                try:
                    ptr = 1
                    # Q values
                    for c in range(nchains):
                        Qchain[station_idx, c] = float(parts[ptr + c]) if ptr + c < len(parts) else 0.0
                    ptr += nchains + 1  # skip aggregate

                    # R values -> W
                    for c in range(nchains):
                        Wchain[station_idx, c] = float(parts[ptr + c]) if ptr + c < len(parts) else 0.0
                    ptr += nchains + 1

                    # U values
                    for c in range(nchains):
                        Uchain[station_idx, c] = float(parts[ptr + c]) if ptr + c < len(parts) else 0.0
                    ptr += nchains + 1

                    # X values -> T
                    for c in range(nchains):
                        Tchain[station_idx, c] = float(parts[ptr + c]) if ptr + c < len(parts) else 0.0
                except (ValueError, IndexError):
                    pass

    return Uchain, Qchain, Wchain, Tchain


def solver_qns(
    sn: NetworkStruct,
    options: Optional[SolverQNSOptions] = None
) -> SolverQNSReturn:
    """
    QNS solver handler - calls external qnsolver via subprocess.

    Performs queueing network analysis using the external qnsolver tool by:
    1. Writing the model to JMVA XML format
    2. Calling qnsolver via subprocess
    3. Parsing the results

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverQNSReturn with all performance metrics

    Raises:
        RuntimeError: If qnsolver is not available or fails
    """
    start_time = time.time()

    if options is None:
        options = SolverQNSOptions()

    if not is_qns_available():
        raise RuntimeError(
            "SolverQNS requires the 'qnsolver' command to be available in your system PATH.\n"
            "You can install it from: http://www.sce.carleton.ca/rads/lqns/"
        )

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    actual_method = options.method

    # Map method to multiserver configuration
    method_map = {
        'conway': 'conway',
        'rolia': 'rolia',
        'zhou': 'zhou',
        'suri': 'suri',
        'reiser': 'reiser',
        'schmidt': 'schmidt',
    }
    if options.method in method_map:
        options.multiserver = method_map[options.method]

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix='qns_')

    try:
        model_path = os.path.join(temp_dir, 'model.jmva')
        result_path = os.path.join(temp_dir, 'result.jmva')
        log_path = os.path.join(temp_dir, 'console.out')

        # Write model to JMVA format
        _write_jmva_file(sn, model_path, options)

        # Build command
        cmd = ['qnsolver', '-l', model_path]

        # Add multiserver method if needed
        if _has_multi_server(sn) and options.multiserver:
            method_flags = {
                'default': '-mconway',
                'conway': '-mconway',
                'reiser': '-mreiser',
                'rolia': '-mrolia',
                'zhou': '-mzhou',
            }
            if options.multiserver in method_flags:
                cmd.append(method_flags[options.multiserver])
                if options.multiserver != 'default':
                    actual_method = options.multiserver
                else:
                    actual_method = 'conway'

        cmd.extend(['-o', result_path])

        if options.verbose:
            print(f"SolverQNS command: {' '.join(cmd)}")

        # Execute command
        with open(log_path, 'w') as log_file:
            result = subprocess.run(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=temp_dir,
                timeout=300
            )

        if result.returncode != 0:
            log_content = ''
            if os.path.exists(log_path):
                with open(log_path, 'r') as f:
                    log_content = f.read()
            raise RuntimeError(f"QNS solver failed with exit code: {result.returncode}\nLog: {log_content}")

        # Parse results
        Uchain, Qchain, Wchain, Tchain = _parse_results(result_path, sn, C)

        # Get chain demands for deaggregation
        demands = sn_get_demands_chain(sn)

        # Calculate system throughput for each chain
        Xchain = np.zeros((1, C))
        refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(K, dtype=int)
        for c in range(C):
            inchain = sn.inchain.get(c, np.array([]))
            if len(inchain) > 0:
                ref_class_idx = int(inchain[0])
                if ref_class_idx < len(refstat):
                    ref_stat_idx = int(refstat[ref_class_idx])
                    Xchain[0, c] = Tchain[ref_stat_idx, c]

        # Adjust utilizations for multi-server stations
        if sn.nservers is not None:
            for i in range(M):
                servers = sn.nservers[i] if len(sn.nservers.shape) == 1 else sn.nservers[i, 0]
                if not np.isinf(servers) and servers > 1:
                    for c in range(C):
                        Uchain[i, c] = Uchain[i, c] / servers

        # Deaggregate chain results to class results
        Rchain = np.nan_to_num(Wchain, nan=0.0)

        deagg = sn_deaggregate_chain_results(
            sn,
            demands.Lchain,
            None,
            demands.STchain,
            demands.Vchain,
            demands.alpha,
            Qchain,
            Uchain,
            Rchain,
            Tchain,
            None,
            Xchain
        )

        runtime = time.time() - start_time

        return SolverQNSReturn(
            Q=deagg.Q,
            U=deagg.U,
            R=deagg.R,
            T=deagg.T,
            A=np.zeros((M, K)),  # Will be calculated later
            W=deagg.R,
            C=deagg.C,
            X=deagg.X,
            runtime=runtime,
            method=actual_method,
            it=0
        )

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
