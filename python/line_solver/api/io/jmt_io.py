"""
JMT (Java Modelling Tools) I/O functions.

This module provides functions for importing from and exporting to JMT file formats:
- JSIMG: JMT simulation model format (XML)
- JMVA: JMT MVA model format (XML)
- JSIM: JMT simulation format

Port from:
    - matlab/src/io/JMT2LINE.m
    - matlab/src/io/JMVA2LINE.m
    - matlab/src/io/JSIM2LINE.m
    - matlab/src/io/QN2JSIMG.m
"""

import os
import tempfile
import xml.etree.ElementTree as ET
from xml.dom import minidom
from typing import Any, Optional, Dict, List, Tuple, Union
from dataclasses import dataclass, field
import numpy as np

from .logging import line_warning, line_error


@dataclass
class JSIMGInfo:
    """Information extracted from a JSIMG file."""
    name: str
    nodes: List[Dict[str, Any]]
    classes: List[Dict[str, Any]]
    connections: List[Tuple[str, str]]
    parameters: Dict[str, Any]


def jmt2line(filename: str) -> Dict[str, Any]:
    """
    Import a JMT model file and convert to LINE network specification.

    Handles JSIMG (simulation) and JMVA (analytical) formats.

    Args:
        filename: Path to JMT file (.jsimg or .jmva)

    Returns:
        Dictionary with network model specification

    References:
        MATLAB: matlab/src/io/JMT2LINE.m
    """
    ext = os.path.splitext(filename)[1].lower()

    if ext == '.jsimg':
        return jsim2line(filename)
    elif ext == '.jmva':
        return jmva2line(filename)
    else:
        # Try to detect format from file content
        try:
            tree = ET.parse(filename)
            root = tree.getroot()
            if root.tag == 'sim' or 'simType' in root.attrib:
                return jsim2line(filename)
            elif root.tag == 'model' or 'modelType' in root.attrib:
                return jmva2line(filename)
            else:
                line_warning('jmt2line', f'Unknown JMT format in {filename}, trying JSIM')
                return jsim2line(filename)
        except Exception as e:
            line_error('jmt2line', f'Failed to parse JMT file: {e}')
            return {}


def jsim2line(filename: str) -> Dict[str, Any]:
    """
    Import a JSIM model file and convert to LINE network specification.

    JSIM files are XML-based simulation model definitions used by JMT.

    Args:
        filename: Path to JSIM file

    Returns:
        Dictionary with network model specification:
            - name: Model name
            - nodes: List of node specifications
            - classes: List of class specifications
            - routing: Routing matrix
            - parameters: Additional parameters

    References:
        MATLAB: matlab/src/io/JSIM2LINE.m
    """
    result = {
        'name': 'imported_model',
        'nodes': [],
        'classes': [],
        'routing': {},
        'parameters': {},
    }

    try:
        tree = ET.parse(filename)
        root = tree.getroot()

        # Get model name
        if 'name' in root.attrib:
            result['name'] = root.attrib['name']
        elif root.find('.//parameters') is not None:
            params = root.find('.//parameters')
            if 'name' in params.attrib:
                result['name'] = params.attrib['name']

        # Parse nodes
        nodes_elem = root.find('.//nodes') or root.findall('.//node')
        if nodes_elem is not None:
            for node_elem in nodes_elem if isinstance(nodes_elem, list) else nodes_elem:
                node_spec = _parse_jsim_node(node_elem)
                if node_spec:
                    result['nodes'].append(node_spec)

        # Parse classes
        classes_elem = root.find('.//userClasses') or root.findall('.//userClass')
        if classes_elem is not None:
            for class_elem in classes_elem if isinstance(classes_elem, list) else classes_elem:
                class_spec = _parse_jsim_class(class_elem)
                if class_spec:
                    result['classes'].append(class_spec)

        # Invert priorities from JMT convention to LINE convention
        _invert_jmt_priorities(result['classes'])

        # Parse connections/routing
        connections_elem = root.find('.//connections') or root.findall('.//connection')
        if connections_elem is not None:
            for conn_elem in connections_elem if isinstance(connections_elem, list) else connections_elem:
                _parse_jsim_connection(conn_elem, result)

        # Parse measures
        measures_elem = root.find('.//measures')
        if measures_elem is not None:
            result['parameters']['measures'] = _parse_jsim_measures(measures_elem)

    except ET.ParseError as e:
        line_error('jsim2line', f'XML parse error: {e}')
    except Exception as e:
        line_error('jsim2line', f'Error importing JSIM: {e}')

    return result


def jmva2line(filename: str) -> Dict[str, Any]:
    """
    Import a JMVA model file and convert to LINE network specification.

    JMVA files are XML-based analytical model definitions used by JMT.

    Args:
        filename: Path to JMVA file

    Returns:
        Dictionary with network model specification

    References:
        MATLAB: matlab/src/io/JMVA2LINE.m
    """
    result = {
        'name': 'imported_mva_model',
        'nodes': [],
        'classes': [],
        'routing': {},
        'parameters': {},
    }

    try:
        tree = ET.parse(filename)
        root = tree.getroot()

        # Get model name
        if 'name' in root.attrib:
            result['name'] = root.attrib['name']

        # Parse stations
        stations_elem = root.find('.//stations') or root.findall('.//station')
        if stations_elem is not None:
            for station_elem in stations_elem if isinstance(stations_elem, list) else stations_elem:
                station_spec = _parse_jmva_station(station_elem)
                if station_spec:
                    result['nodes'].append(station_spec)

        # Parse classes
        classes_elem = root.find('.//classes') or root.findall('.//class')
        if classes_elem is not None:
            for class_elem in classes_elem if isinstance(classes_elem, list) else classes_elem:
                class_spec = _parse_jmva_class(class_elem)
                if class_spec:
                    result['classes'].append(class_spec)

        # Parse service demands
        demands_elem = root.find('.//serviceDemands') or root.find('.//demands')
        if demands_elem is not None:
            result['parameters']['demands'] = _parse_jmva_demands(demands_elem)

        # Parse visits
        visits_elem = root.find('.//visits')
        if visits_elem is not None:
            result['parameters']['visits'] = _parse_jmva_visits(visits_elem)

    except ET.ParseError as e:
        line_error('jmva2line', f'XML parse error: {e}')
    except Exception as e:
        line_error('jmva2line', f'Error importing JMVA: {e}')

    return result


def qn2jsimg(model: Any, output_filename: Optional[str] = None,
             options: Optional[Dict] = None) -> str:
    """
    Export a Network model to JMT JSIMG format.

    Creates a JSIMG (JMT simulation) XML file from a Network model.

    Args:
        model: Network model or NetworkStruct
        output_filename: Optional output file path (default: temp file)
        options: Optional export options

    Returns:
        Path to the created JSIMG file

    References:
        MATLAB: matlab/src/io/QN2JSIMG.m
    """
    if output_filename is None:
        fd, output_filename = tempfile.mkstemp(suffix='.jsimg')
        os.close(fd)

    if options is None:
        options = {}

    # Get network structure
    if hasattr(model, 'getStruct'):
        sn = model.getStruct()
    else:
        sn = model

    # Build XML
    root = ET.Element('archive')
    root.set('name', getattr(sn, 'name', 'model'))
    root.set('timestamp', '')
    root.set('xsi:noNamespaceSchemaLocation', 'Archive.xsd')
    root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')

    sim = ET.SubElement(root, 'sim')
    sim.set('name', getattr(sn, 'name', 'model'))
    sim.set('xsi:noNamespaceSchemaLocation', 'SIMmodeldefinition.xsd')
    sim.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')

    # Add simulation parameters
    _add_jsimg_parameters(sim, sn, options)

    # Add user classes
    _add_jsimg_classes(sim, sn)

    # Add nodes
    _add_jsimg_nodes(sim, sn)

    # Add measures
    _add_jsimg_measures(sim, sn)

    # Add connections
    _add_jsimg_connections(sim, sn)

    # Write XML file
    xml_str = ET.tostring(root, encoding='unicode')
    # Pretty print
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')

    with open(output_filename, 'w') as f:
        f.write(pretty_xml)

    return output_filename


# Helper functions for JSIM parsing

def _parse_jsim_node(node_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JSIM node element."""
    node_spec = {
        'name': node_elem.get('name', 'Unknown'),
        'type': 'Queue',
    }

    # Determine node type from className or tag
    class_name = node_elem.get('className', '').lower()
    if 'source' in class_name:
        node_spec['type'] = 'Source'
    elif 'sink' in class_name:
        node_spec['type'] = 'Sink'
    elif 'delay' in class_name or 'infinite' in class_name:
        node_spec['type'] = 'Delay'
    elif 'fork' in class_name:
        node_spec['type'] = 'Fork'
    elif 'join' in class_name:
        node_spec['type'] = 'Join'
    elif 'router' in class_name:
        node_spec['type'] = 'Router'
    elif 'queue' in class_name or 'server' in class_name:
        node_spec['type'] = 'Queue'

    # Parse section parameters
    for section in node_elem.findall('.//section'):
        section_class = section.get('className', '').lower()
        if 'queue' in section_class:
            # Parse queue parameters
            for param in section.findall('.//parameter'):
                param_name = param.get('name', '').lower()
                if 'size' in param_name or 'capacity' in param_name:
                    try:
                        value = param.find('.//value')
                        if value is not None and value.text:
                            node_spec['capacity'] = int(value.text)
                    except (ValueError, TypeError):
                        pass
                elif 'server' in param_name:
                    try:
                        value = param.find('.//value')
                        if value is not None and value.text:
                            node_spec['servers'] = int(value.text)
                    except (ValueError, TypeError):
                        pass

        elif 'service' in section_class:
            # Parse service parameters
            services = {}
            for subparam in section.findall('.//subParameter'):
                class_name = subparam.get('classPath', '')
                if class_name:
                    service_info = {}
                    for param in subparam.findall('.//parameter'):
                        param_name = param.get('name', '').lower()
                        value_elem = param.find('.//value')
                        if value_elem is not None and value_elem.text:
                            try:
                                if 'mean' in param_name or 'rate' in param_name:
                                    service_info[param_name] = float(value_elem.text)
                            except ValueError:
                                pass
                    services[class_name] = service_info
            if services:
                node_spec['services'] = services

    return node_spec


def _parse_jsim_class(class_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JSIM user class element."""
    class_spec = {
        'name': class_elem.get('name', 'Unknown'),
        'type': 'closed',
    }

    class_type = class_elem.get('type', '').lower()
    if 'open' in class_type:
        class_spec['type'] = 'open'
        class_spec['population'] = float('inf')
    else:
        class_spec['type'] = 'closed'
        try:
            pop = class_elem.get('population', '1')
            class_spec['population'] = int(pop)
        except ValueError:
            class_spec['population'] = 1

    # Get reference station
    ref_station = class_elem.get('referenceStation', '')
    if ref_station:
        class_spec['refstation'] = ref_station

    # Get priority (will be inverted later in _invert_jmt_priorities)
    # JMT uses higher value = higher priority, LINE uses lower value = higher priority
    try:
        prio = class_elem.get('priority', '0')
        class_spec['priority'] = int(prio)
    except ValueError:
        class_spec['priority'] = 0

    return class_spec


def _invert_jmt_priorities(classes: List[Dict[str, Any]]) -> None:
    """Invert priorities from JMT convention to LINE convention.

    JMT uses higher priority value = higher priority.
    LINE uses lower priority value = higher priority.

    Args:
        classes: List of class specifications with 'priority' field
    """
    if not classes:
        return

    # Find max priority
    max_prio = 0
    for cls in classes:
        prio = cls.get('priority', 0)
        if prio > max_prio:
            max_prio = prio

    # Invert each priority
    for cls in classes:
        jmt_prio = cls.get('priority', 0)
        cls['priority'] = max_prio - jmt_prio


def _parse_jsim_connection(conn_elem: ET.Element, result: Dict) -> None:
    """Parse a JSIM connection element and add to routing."""
    source = conn_elem.get('source', '')
    target = conn_elem.get('target', '')
    if source and target:
        if 'connections' not in result:
            result['connections'] = []
        result['connections'].append((source, target))


def _parse_jsim_measures(measures_elem: ET.Element) -> List[Dict[str, Any]]:
    """Parse JSIM measures elements."""
    measures = []
    for measure in measures_elem.findall('.//measure'):
        measure_spec = {
            'type': measure.get('measureType', ''),
            'station': measure.get('station', ''),
            'class': measure.get('class', ''),
        }
        measures.append(measure_spec)
    return measures


# Helper functions for JMVA parsing

def _parse_jmva_station(station_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JMVA station element."""
    station_spec = {
        'name': station_elem.get('name', 'Unknown'),
        'type': 'Queue',
    }

    station_type = station_elem.get('type', '').lower()
    if 'delay' in station_type or 'li' in station_type:
        station_spec['type'] = 'Delay'
    elif 'ld' in station_type:
        station_spec['type'] = 'Queue'
        station_spec['load_dependent'] = True

    servers = station_elem.get('servers', '1')
    try:
        station_spec['servers'] = int(servers)
    except ValueError:
        station_spec['servers'] = 1

    return station_spec


def _parse_jmva_class(class_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JMVA class element."""
    class_spec = {
        'name': class_elem.get('name', 'Unknown'),
        'type': 'closed',
    }

    class_type = class_elem.get('type', '').lower()
    if 'open' in class_type:
        class_spec['type'] = 'open'
        class_spec['population'] = float('inf')
        rate = class_elem.get('rate', '1.0')
        try:
            class_spec['arrival_rate'] = float(rate)
        except ValueError:
            class_spec['arrival_rate'] = 1.0
    else:
        class_spec['type'] = 'closed'
        pop = class_elem.get('population', '1')
        try:
            class_spec['population'] = int(pop)
        except ValueError:
            class_spec['population'] = 1

    return class_spec


def _parse_jmva_demands(demands_elem: ET.Element) -> Dict[str, Dict[str, float]]:
    """Parse JMVA service demands."""
    demands = {}
    for demand in demands_elem.findall('.//serviceDemand') or demands_elem.findall('.//demand'):
        station = demand.get('stationName', demand.get('station', ''))
        job_class = demand.get('customerClass', demand.get('class', ''))
        value = demand.text or demand.get('value', '0')
        try:
            demand_value = float(value)
            if station not in demands:
                demands[station] = {}
            demands[station][job_class] = demand_value
        except ValueError:
            pass
    return demands


def _parse_jmva_visits(visits_elem: ET.Element) -> Dict[str, Dict[str, float]]:
    """Parse JMVA visits."""
    visits = {}
    for visit in visits_elem.findall('.//visit'):
        station = visit.get('stationName', visit.get('station', ''))
        job_class = visit.get('customerClass', visit.get('class', ''))
        value = visit.text or visit.get('value', '1')
        try:
            visit_value = float(value)
            if station not in visits:
                visits[station] = {}
            visits[station][job_class] = visit_value
        except ValueError:
            pass
    return visits


# Helper functions for JSIMG export

def _add_jsimg_parameters(sim: ET.Element, sn: Any, options: Dict) -> None:
    """Add simulation parameters to JSIMG."""
    params = ET.SubElement(sim, 'parameters')

    # Seed
    seed = ET.SubElement(params, 'seed')
    seed.text = str(options.get('seed', 1))

    # Simulation seconds
    sim_seconds = ET.SubElement(params, 'maxTime')
    sim_seconds.text = str(options.get('max_time', -1))

    # Max samples
    max_samples = ET.SubElement(params, 'maxSamples')
    max_samples.text = str(options.get('max_samples', 1000000))


def _add_jsimg_classes(sim: ET.Element, sn: Any) -> None:
    """Add user classes to JSIMG."""
    # JMT uses higher priority value = higher priority, LINE uses lower value = higher priority
    # We need to invert priorities when exporting to JMT
    max_prio = 0
    if hasattr(sn, 'classprio') and sn.classprio is not None:
        max_prio = int(np.max(sn.classprio))

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'
        njobs = sn.njobs[k] if hasattr(sn, 'njobs') else 0

        user_class = ET.SubElement(sim, 'userClass')
        user_class.set('name', class_name)

        # Set priority with inversion: LINE uses lower=higher, JMT uses higher=higher
        if hasattr(sn, 'classprio') and sn.classprio is not None:
            line_prio = int(sn.classprio[k]) if sn.classprio.ndim == 1 else int(sn.classprio[0, k])
            jmt_prio = max_prio - line_prio
            user_class.set('priority', str(jmt_prio))
        else:
            user_class.set('priority', '0')

        if np.isinf(njobs):
            user_class.set('type', 'open')
        else:
            user_class.set('type', 'closed')
            user_class.set('population', str(int(njobs)))

            # Reference station
            if hasattr(sn, 'refstat') and hasattr(sn, 'stationToNode'):
                refstat = sn.refstat[k]
                ref_node = sn.stationToNode[refstat]
                ref_name = sn.nodenames[ref_node] if hasattr(sn, 'nodenames') else f'Station{refstat}'
                user_class.set('referenceStation', ref_name)


def _add_jsimg_nodes(sim: ET.Element, sn: Any) -> None:
    """Add nodes to JSIMG."""
    for i in range(sn.nnodes):
        node_name = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
        node_type = sn.nodetype[i] if hasattr(sn, 'nodetype') else None
        type_name = node_type.name if hasattr(node_type, 'name') else 'QUEUE'

        node = ET.SubElement(sim, 'node')
        node.set('name', node_name)

        if type_name == 'SOURCE':
            node.set('className', 'RandomSource')
            _add_source_sections(node, sn, i)
        elif type_name == 'SINK':
            node.set('className', 'Sink')
        elif type_name == 'DELAY':
            node.set('className', 'Delay')
            _add_delay_sections(node, sn, i)
        elif type_name == 'QUEUE':
            node.set('className', 'Server')
            _add_queue_sections(node, sn, i)
        elif type_name == 'FORK':
            node.set('className', 'Fork')
        elif type_name == 'JOIN':
            node.set('className', 'Join')
        elif type_name == 'ROUTER':
            node.set('className', 'Router')


def _add_source_sections(node: ET.Element, sn: Any, node_idx: int) -> None:
    """Add source node sections."""
    section = ET.SubElement(node, 'section')
    section.set('className', 'RandomSource')

    ist = sn.nodeToStation[node_idx] if hasattr(sn, 'nodeToStation') else node_idx

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        # Get arrival rate
        rate = 0
        if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1]:
            rate = sn.rates[ist, k]
            if np.isnan(rate):
                rate = 0

        if rate > 0:
            param = ET.SubElement(section, 'parameter')
            param.set('classPath', class_name)
            param.set('name', 'ServiceStrategy')

            subparam = ET.SubElement(param, 'subParameter')
            subparam.set('classPath', 'jmt.engine.random.Exponential')
            subparam.set('name', 'Exponential')

            rate_param = ET.SubElement(subparam, 'parameter')
            rate_param.set('name', 'lambda')
            value = ET.SubElement(rate_param, 'value')
            value.text = str(rate)


def _add_delay_sections(node: ET.Element, sn: Any, node_idx: int) -> None:
    """Add delay node sections."""
    section = ET.SubElement(node, 'section')
    section.set('className', 'ServiceSection')

    ist = sn.nodeToStation[node_idx] if hasattr(sn, 'nodeToStation') else node_idx

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        # Get service rate
        rate = 0
        if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1]:
            rate = sn.rates[ist, k]
            if np.isnan(rate):
                continue

        mean = 1.0 / rate if rate > 0 else 0

        param = ET.SubElement(section, 'parameter')
        param.set('classPath', class_name)
        param.set('name', 'ServiceStrategy')

        subparam = ET.SubElement(param, 'subParameter')
        subparam.set('classPath', 'jmt.engine.random.Exponential')
        subparam.set('name', 'Exponential')

        mean_param = ET.SubElement(subparam, 'parameter')
        mean_param.set('name', 'mean')
        value = ET.SubElement(mean_param, 'value')
        value.text = str(mean)


def _add_queue_sections(node: ET.Element, sn: Any, node_idx: int) -> None:
    """Add queue node sections."""
    ist = sn.nodeToStation[node_idx] if hasattr(sn, 'nodeToStation') else node_idx

    # Queue section
    queue_section = ET.SubElement(node, 'section')
    queue_section.set('className', 'Queue')

    # Server count
    servers = 1
    if hasattr(sn, 'nservers') and ist < len(sn.nservers):
        servers = int(sn.nservers[ist])

    servers_param = ET.SubElement(queue_section, 'parameter')
    servers_param.set('name', 'size')
    value = ET.SubElement(servers_param, 'value')
    value.text = str(servers)

    # Service section
    service_section = ET.SubElement(node, 'section')
    service_section.set('className', 'ServiceSection')

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        # Get service rate
        rate = 0
        if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1]:
            rate = sn.rates[ist, k]
            if np.isnan(rate):
                continue

        mean = 1.0 / rate if rate > 0 else 0

        param = ET.SubElement(service_section, 'parameter')
        param.set('classPath', class_name)
        param.set('name', 'ServiceStrategy')

        subparam = ET.SubElement(param, 'subParameter')
        subparam.set('classPath', 'jmt.engine.random.Exponential')
        subparam.set('name', 'Exponential')

        mean_param = ET.SubElement(subparam, 'parameter')
        mean_param.set('name', 'mean')
        value = ET.SubElement(mean_param, 'value')
        value.text = str(mean)


def _add_jsimg_measures(sim: ET.Element, sn: Any) -> None:
    """Add measures to JSIMG."""
    measures = ET.SubElement(sim, 'measures')

    for ist in range(sn.nstations):
        if hasattr(sn, 'stationToNode'):
            node_idx = sn.stationToNode[ist]
        else:
            node_idx = ist
        station_name = sn.nodenames[node_idx] if hasattr(sn, 'nodenames') else f'Station{ist}'

        for k in range(sn.nclasses):
            class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

            # Queue length
            measure = ET.SubElement(measures, 'measure')
            measure.set('measureType', 'Queue length')
            measure.set('station', station_name)
            measure.set('class', class_name)

            # Throughput
            measure = ET.SubElement(measures, 'measure')
            measure.set('measureType', 'Throughput')
            measure.set('station', station_name)
            measure.set('class', class_name)


def _add_jsimg_connections(sim: ET.Element, sn: Any) -> None:
    """Add connections to JSIMG."""
    if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
        for k in range(sn.nclasses):
            for i in range(sn.nnodes):
                for j in range(sn.nnodes):
                    idx_from = i * sn.nclasses + k
                    idx_to = j * sn.nclasses + k
                    if idx_from < sn.rtnodes.shape[0] and idx_to < sn.rtnodes.shape[1]:
                        prob = sn.rtnodes[idx_from, idx_to]
                        if prob > 0:
                            conn = ET.SubElement(sim, 'connection')
                            source_name = sn.nodenames[i] if hasattr(sn, 'nodenames') else f'Node{i}'
                            target_name = sn.nodenames[j] if hasattr(sn, 'nodenames') else f'Node{j}'
                            conn.set('source', source_name)
                            conn.set('target', target_name)


__all__ = [
    'JSIMGInfo',
    'jmt2line',
    'jsim2line',
    'jmva2line',
    'qn2jsimg',
]
