"""
JMVA XML file writer for qnsolver.

This module generates JMVA (JMT MVA) format XML files from NetworkStruct,
which can be processed by the qnsolver command-line tool.
"""

import numpy as np
from xml.etree.ElementTree import Element, SubElement, ElementTree
from xml.dom import minidom
from typing import Dict, Any, Optional
import os

from ...api.sn import NetworkStruct, NodeType
from ...api.sn.demands import sn_get_demands_chain


def write_jmva(sn: NetworkStruct, output_path: str, options: Optional[Dict[str, Any]] = None) -> str:
    """
    Write a NetworkStruct to JMVA XML format.

    Args:
        sn: NetworkStruct containing the queueing network model
        output_path: Path to write the JMVA XML file
        options: Optional solver options (method, samples, etc.)

    Returns:
        Path to the written file
    """
    if options is None:
        options = {}

    method = options.get('method', 'default')
    samples = options.get('samples', 10000)

    # Create root element
    model = Element('model')
    model.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    model.set('xsi:noNamespaceSchemaLocation', 'JMTmodel.xsd')

    # Get demands per chain using proper chain aggregation
    demands = sn_get_demands_chain(sn)
    Lchain = demands.Lchain      # Dchain in JMVA context = Lchain
    STchain = demands.STchain
    Vchain = demands.Vchain
    Nchain = demands.Nchain.flatten()

    # Parameters element
    parameters = SubElement(model, 'parameters')

    # Classes element
    classes = SubElement(parameters, 'classes')
    classes.set('number', str(sn.nchains))

    # Determine source nodes
    sourceid = []
    for i, nt in enumerate(sn.nodetype):
        sourceid.append(nt == NodeType.SOURCE)

    # Add class elements
    for c in range(sn.nchains):
        # Sum of jobs in this chain
        sum_njobs = 0.0
        if sn.chains is not None and c < sn.chains.shape[0]:
            for k in range(sn.chains.shape[1]):
                if sn.chains[c, k] > 0 and k < len(sn.njobs):
                    sum_njobs += sn.njobs[k]

        if np.isfinite(sum_njobs) and not np.isnan(sum_njobs):
            # Closed class
            class_elem = SubElement(classes, 'closedclass')
            class_elem.set('population', str(int(Nchain[c])))
            class_elem.set('name', f'Chain{c+1:02d}')
        else:
            # Open class - calculate arrival rate
            rate_sum = 0.0
            for i, is_source in enumerate(sourceid):
                if is_source:
                    if sn.chains is not None and c < sn.chains.shape[0]:
                        for k in range(sn.chains.shape[1]):
                            if sn.chains[c, k] > 0 and k < sn.rates.shape[1]:
                                rate_sum += sn.rates[i, k]
            class_elem = SubElement(classes, 'openclass')
            class_elem.set('rate', str(rate_sum))
            class_elem.set('name', f'Chain{c+1:02d}')

    # Count stations (excluding sources)
    num_stations = sn.nstations
    for nt in sn.nodetype:
        if nt == NodeType.SOURCE:
            num_stations -= 1

    # Stations element
    stations = SubElement(parameters, 'stations')
    stations.set('number', str(num_stations))

    # Track load-dependent stations
    is_load_dep = [False] * sn.nstations

    # Add station elements
    for i in range(sn.nstations):
        node_idx = int(sn.stationToNode[i])
        node_type = sn.nodetype[node_idx]

        if node_type == NodeType.DELAY:
            stat_elem = SubElement(stations, 'delaystation')
            stat_elem.set('name', sn.nodenames[node_idx])
        elif node_type == NodeType.QUEUE:
            nservers = sn.nservers[i] if i < len(sn.nservers) else 1
            if nservers == 1:
                is_load_dep[i] = False
                stat_elem = SubElement(stations, 'listation')
            else:
                is_load_dep[i] = True
                stat_elem = SubElement(stations, 'ldstation')
            stat_elem.set('name', sn.nodenames[node_idx])
            stat_elem.set('servers', '1')
        elif node_type == NodeType.SOURCE:
            continue  # Skip sources
        else:
            continue  # Skip other types

        # Service times
        srv_times = SubElement(stat_elem, 'servicetimes')
        for c in range(sn.nchains):
            st_value = STchain[i, c] if i < STchain.shape[0] and c < STchain.shape[1] else 0.0
            if is_load_dep[i]:
                # Load-dependent station
                stat_srv_time = SubElement(srv_times, 'servicetimes')
                stat_srv_time.set('customerclass', f'Chain{c+1:02d}')

                # Build load-dependent service time string
                # Position k = service time when k customers present
                # For multiserver: ST(n) = ST / min(n, nservers)
                total_pop = int(np.sum([n for n in Nchain if np.isfinite(n)]))
                nservers = sn.nservers[i] if i < len(sn.nservers) else 1
                ld_srv_string = str(st_value)  # n=1: base service time
                for n in range(2, total_pop + 1):
                    ld_srv_string += f';{st_value / min(n, nservers)}'
                stat_srv_time.text = ld_srv_string
            else:
                stat_srv_time = SubElement(srv_times, 'servicetime')
                stat_srv_time.set('customerclass', f'Chain{c+1:02d}')
                stat_srv_time.text = str(st_value)

        # Visits
        visits = SubElement(stat_elem, 'visits')
        for c in range(sn.nchains):
            visit_elem = SubElement(visits, 'visit')
            visit_elem.set('customerclass', f'Chain{c+1:02d}')

            st_value = STchain[i, c] if i < STchain.shape[0] and c < STchain.shape[1] else 0.0
            l_value = Lchain[i, c] if i < Lchain.shape[0] and c < Lchain.shape[1] else 0.0

            if st_value > 0:
                val = l_value / st_value
            else:
                val = 0.0
            visit_elem.text = str(val)

    # Reference stations
    ref_stations = SubElement(parameters, 'ReferenceStation')
    ref_stations.set('number', str(sn.nchains))

    for c in range(sn.nchains):
        class_ref = SubElement(ref_stations, 'Class')
        class_ref.set('name', f'Chain{c+1:02d}')

        # Get reference station for this chain
        if sn.inchain is not None and c < len(sn.inchain):
            inchain = sn.inchain[c]
            if len(inchain) > 0:
                first_class = int(inchain[0])
                if sn.refstat is not None and first_class < len(sn.refstat):
                    refstat_idx = int(sn.refstat[first_class])
                    node_idx = int(sn.stationToNode[refstat_idx])
                    class_ref.set('refStation', sn.nodenames[node_idx])
                else:
                    class_ref.set('refStation', sn.nodenames[int(sn.stationToNode[0])])
            else:
                class_ref.set('refStation', sn.nodenames[int(sn.stationToNode[0])])
        else:
            class_ref.set('refStation', sn.nodenames[int(sn.stationToNode[0])])

    # Algorithm parameters
    alg_params = SubElement(model, 'algParams')
    alg_type = SubElement(alg_params, 'algType')
    _set_algorithm_name(alg_type, sn, method)
    alg_type.set('tolerance', '1.0E-7')
    alg_type.set('maxSamples', str(samples))
    compare_algs = SubElement(alg_params, 'compareAlgs')
    compare_algs.set('value', 'false')

    # Write to file
    tree = ElementTree(model)

    # Pretty print
    xml_str = minidom.parseString(
        _element_to_string(model)
    ).toprettyxml(indent='  ')

    # Remove extra blank lines
    lines = [line for line in xml_str.split('\n') if line.strip()]
    xml_str = '\n'.join(lines)

    with open(output_path, 'w') as f:
        f.write(xml_str)

    return output_path


def _set_algorithm_name(alg_elem: Element, sn: NetworkStruct, method: str) -> None:
    """Set algorithm name based on method and network characteristics."""
    method = method.lower()

    # Check for multi-server
    has_multi_server = False
    if sn.nservers is not None:
        for ns in sn.nservers:
            if ns > 1 and not np.isinf(ns):
                has_multi_server = True
                break

    if method == 'jmva.recal':
        alg_elem.set('name', 'RECAL')
    elif method == 'jmva.comom':
        alg_elem.set('name', 'CoMoM')
    elif method == 'jmva.chow':
        alg_elem.set('name', 'Chow')
    elif method in ('jmva.bs', 'jmva.amva'):
        alg_elem.set('name', 'Bard-Schweitzer')
    elif method == 'jmva.aql':
        alg_elem.set('name', 'AQL')
    elif method == 'jmva.lin':
        alg_elem.set('name', 'Linearizer')
    elif method == 'jmva.dmlin':
        alg_elem.set('name', 'De Souza-Muntz Linearizer')
    else:
        alg_elem.set('name', 'MVA')


def _element_to_string(elem: Element) -> bytes:
    """Convert Element to string without XML declaration."""
    from xml.etree.ElementTree import tostring
    return tostring(elem, encoding='unicode').encode('utf-8')
