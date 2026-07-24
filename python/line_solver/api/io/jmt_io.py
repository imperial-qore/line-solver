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

        # The authoritative model lives under <sim>; the sibling <jmodel> block
        # only carries GUI data (station positions, userClass colors) and MUST
        # NOT be parsed -- its duplicate <userClass name=...> entries would
        # otherwise double every class. Scope all lookups to <sim>.
        sim = root.find('sim')
        if sim is None:
            sim = root.find('.//sim')
        if sim is None:
            sim = root

        # Get model name
        if 'name' in sim.attrib:
            result['name'] = sim.attrib['name']
        elif 'name' in root.attrib:
            result['name'] = root.attrib['name']

        # Parse nodes (direct <node> children of <sim>)
        for node_elem in sim.findall('node'):
            node_spec = _parse_jsim_node(node_elem)
            if node_spec:
                result['nodes'].append(node_spec)

        # Parse classes (direct <userClass> children of <sim>)
        for class_elem in sim.findall('userClass'):
            class_spec = _parse_jsim_class(class_elem)
            if class_spec:
                result['classes'].append(class_spec)

        # Invert priorities from JMT convention to LINE convention
        _invert_jmt_priorities(result['classes'])

        # Parse connections/routing (direct <connection> children of <sim>)
        for conn_elem in sim.findall('connection'):
            _parse_jsim_connection(conn_elem, result)

        # Parse measures
        measures_elem = sim.find('measures')
        if measures_elem is not None:
            result['parameters']['measures'] = _parse_jsim_measures(measures_elem)

        # Parse the initial SPN marking (<preload>): the per-place, per-class
        # token counts used to seed the initial state (mirrors the Preload
        # branch of matlab/src/io/JSIM2LINE.m).
        preload_elem = sim.find('preload')
        if preload_elem is not None:
            preload: Dict[str, Dict[str, int]] = {}
            for sp in preload_elem.findall('stationPopulations'):
                station = sp.get('stationName', '')
                if not station:
                    continue
                per_class: Dict[str, int] = {}
                for cp in sp.findall('classPopulation'):
                    cls = cp.get('refClass', '')
                    try:
                        per_class[cls] = int(cp.get('population', '0'))
                    except ValueError:
                        per_class[cls] = 0
                preload[station] = per_class
            if preload:
                result['preload'] = preload

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

    # Add finite capacity regions (blocking regions)
    _add_jsimg_regions(sim, model, sn)

    # Write XML file
    xml_str = ET.tostring(root, encoding='unicode')
    # Pretty print
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')

    with open(output_filename, 'w') as f:
        f.write(pretty_xml)

    return output_filename


# Helper functions for JSIM parsing

def _walk_class_array(param_elem: ET.Element):
    """Yield (className, subParameter) pairs from a JMT per-class array parameter.

    JMT lays out per-class array parameters (ServiceStrategy, RoutingStrategy,
    ...) as an interleaved sequence of <refClass> tags and the <subParameter>
    that applies to the class named by the preceding <refClass>.
    """
    current = None
    for child in list(param_elem):
        if child.tag == 'refClass':
            current = (child.text or '').strip()
        elif child.tag == 'subParameter' and current is not None:
            yield current, child
            current = None


def _find_param(section: ET.Element, name: str) -> Optional[ET.Element]:
    """Return the direct <parameter name=...> child of a section, if present."""
    for param in section.findall('parameter'):
        if param.get('name') == name:
            return param
    return None


def _parse_jmt_distr(strat_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JMT ServiceTimeStrategy subParameter into a distribution spec.

    Returns {'type': <jmt distr name>, 'params': {name: float}} or None when the
    strategy is empty (encoded by JMT as <value>null</value> with no nested
    distribution), which denotes a disabled arrival/service for that class.
    """
    subs = strat_elem.findall('subParameter')
    if not subs:
        return None
    dtype = subs[0].get('name', '')
    params = {}
    if len(subs) > 1:
        # The second subParameter is distrPar, holding the named scalar values.
        for p in subs[1].findall('.//subParameter'):
            pname = p.get('name', '')
            value_elem = p.find('value')
            if pname and value_elem is not None and value_elem.text:
                try:
                    params[pname] = float(value_elem.text)
                except ValueError:
                    pass
    return {'type': dtype, 'params': params}


def _parse_service_array(section: ET.Element) -> Dict[str, Optional[Dict[str, Any]]]:
    """Parse the per-class ServiceStrategy array of a RandomSource/Server/Delay
    section into {className: distr_spec_or_None}."""
    result = {}
    param = _find_param(section, 'ServiceStrategy')
    if param is None:
        return result
    for class_name, strat in _walk_class_array(param):
        result[class_name] = _parse_jmt_distr(strat)
    return result


def _parse_routing_array(section: ET.Element) -> Dict[str, Dict[str, Any]]:
    """Parse the per-class RoutingStrategy array of a Router section into
    {className: {'strategy': <jmt name>, 'dests': {stationName: prob}}}."""
    result = {}
    param = _find_param(section, 'RoutingStrategy')
    if param is None:
        return result
    for class_name, strat in _walk_class_array(param):
        entry = {'strategy': strat.get('name', 'Random'), 'dests': {}}
        for emp in strat.findall('.//subParameter[@name="EmpiricalEntry"]'):
            st = emp.find('.//subParameter[@name="stationName"]/value')
            pr = emp.find('.//subParameter[@name="probability"]/value')
            if st is not None and st.text and pr is not None and pr.text:
                try:
                    entry['dests'][st.text.strip()] = float(pr.text)
                except ValueError:
                    pass
        result[class_name] = entry
    return result


def _parse_transition_vectors(matrix_elem: ET.Element) -> Dict[str, Dict[str, int]]:
    """Parse one TransitionMatrix (a single mode's enabling/inhibiting/firing
    condition) into {stationName: {className: tokens}}.

    A TransitionMatrix wraps one array subParameter of TransitionVectors; each
    vector carries a stationName value plus a per-class array of integer entries
    (interleaved refClass/subParameter, as elsewhere in JMT). Mirrors the
    enablingVector/firingVector traversal in matlab/src/io/JSIM2LINE.m.
    """
    result: Dict[str, Dict[str, int]] = {}
    # The matrix has exactly one child array subParameter holding the vectors.
    vectors_arr = None
    for sub in matrix_elem.findall('subParameter'):
        vectors_arr = sub
        break
    if vectors_arr is None:
        return result
    for vec in vectors_arr.findall('subParameter'):
        station = None
        entries: Dict[str, int] = {}
        for s in vec.findall('subParameter'):
            if s.get('name') == 'stationName':
                v = s.find('value')
                station = v.text.strip() if v is not None and v.text else None
            else:
                # Per-class entries array (enablingEntries/inhibitingEntries/
                # firingEntries): interleaved refClass + entry subParameters.
                for cls, entry in _walk_class_array(s):
                    v = entry.find('value')
                    if v is not None and v.text:
                        try:
                            entries[cls] = int(v.text)
                        except ValueError:
                            pass
        if station is not None:
            result[station] = entries
    return result


def _mode_matrices(section: Optional[ET.Element], param_name: str) -> List[Dict[str, Dict[str, int]]]:
    """Parse a per-mode array of TransitionMatrix parameters (enablingConditions,
    inhibitingConditions, firingOutcomes) into a list indexed by mode."""
    result: List[Dict[str, Dict[str, int]]] = []
    if section is None:
        return result
    param = _find_param(section, param_name)
    if param is None:
        return result
    for mode_matrix in param.findall('subParameter'):
        result.append(_parse_transition_vectors(mode_matrix))
    return result


def _mode_scalars(section: Optional[ET.Element], param_name: str) -> List[Optional[float]]:
    """Parse a per-mode array of scalar-valued parameters (numbersOfServers,
    firingPriorities, firingWeights) into a list indexed by mode."""
    result: List[Optional[float]] = []
    if section is None:
        return result
    param = _find_param(section, param_name)
    if param is None:
        return result
    for sub in param.findall('subParameter'):
        v = sub.find('value')
        if v is not None and v.text:
            txt = v.text.strip()
            try:
                result.append(int(txt))
            except ValueError:
                try:
                    result.append(float(txt))
                except ValueError:
                    result.append(None)
        else:
            result.append(None)
    return result


def _mode_timing(section: Optional[ET.Element]) -> List[Tuple[str, Optional[Dict[str, Any]]]]:
    """Parse the per-mode timingStrategies array into a list of
    (kind, distr_spec) where kind is 'IMMEDIATE' (ZeroServiceTimeStrategy, spec
    None) or 'TIMED' (spec from _parse_jmt_distr). Mirrors the timing branch of
    matlab/src/io/JSIM2LINE.m."""
    result: List[Tuple[str, Optional[Dict[str, Any]]]] = []
    if section is None:
        return result
    param = _find_param(section, 'timingStrategies')
    if param is None:
        return result
    for sub in param.findall('subParameter'):
        cp = sub.get('classPath', '')
        if cp.endswith('ZeroServiceTimeStrategy'):
            result.append(('IMMEDIATE', None))
        else:
            result.append(('TIMED', _parse_jmt_distr(sub)))
    return result


def _parse_transition_modes(sections_map: Dict[str, ET.Element]) -> List[Dict[str, Any]]:
    """Combine the Enabling/Timing/Firing sections of a JMT Transition node into
    a per-mode list of firing specifications. Mirrors the Transition branch of
    matlab/src/io/JSIM2LINE.m."""
    enabling_sec = sections_map.get('Enabling')
    timing_sec = sections_map.get('Timing')
    firing_sec = sections_map.get('Firing')

    mode_names: List[str] = []
    if timing_sec is not None:
        mn = _find_param(timing_sec, 'modeNames')
        if mn is not None:
            for sub in mn.findall('subParameter'):
                v = sub.find('value')
                mode_names.append(v.text.strip() if v is not None and v.text else '')
    nmodes = len(mode_names)

    enabling_conds = _mode_matrices(enabling_sec, 'enablingConditions')
    inhibiting_conds = _mode_matrices(enabling_sec, 'inhibitingConditions')
    firing_outcomes = _mode_matrices(firing_sec, 'firingOutcomes')
    servers = _mode_scalars(timing_sec, 'numbersOfServers')
    priorities = _mode_scalars(timing_sec, 'firingPriorities')
    weights = _mode_scalars(timing_sec, 'firingWeights')
    timing_strats = _mode_timing(timing_sec)

    modes: List[Dict[str, Any]] = []
    for m in range(nmodes):
        modes.append({
            'name': mode_names[m],
            'servers': servers[m] if m < len(servers) else -1,
            'firing_priority': priorities[m] if m < len(priorities) else -1,
            'firing_weight': weights[m] if m < len(weights) else 1.0,
            'timing': timing_strats[m] if m < len(timing_strats) else ('TIMED', None),
            'enabling': enabling_conds[m] if m < len(enabling_conds) else {},
            'inhibiting': inhibiting_conds[m] if m < len(inhibiting_conds) else {},
            'firing': firing_outcomes[m] if m < len(firing_outcomes) else {},
        })
    return modes


def _parse_jsim_node(node_elem: ET.Element) -> Optional[Dict[str, Any]]:
    """Parse a JSIM node element.

    Mirrors matlab/src/io/JSIM2LINE.m: the node type is inferred from its section
    classNames (RandomSource -> Source, JobSink -> Sink, Delay -> Delay,
    Server/PSServer -> Queue, Storage -> Place, Enabling/Timing/Firing ->
    Transition, ...), arrival distributions come from the RandomSource section,
    service distributions from the Server/PSServer/Delay section, and per-class
    routing from the Router section.
    """
    node_spec = {
        'name': node_elem.get('name', 'Unknown'),
        'type': 'Queue',
    }

    sections = node_elem.findall('section')
    sec_classes = [s.get('className', '') for s in sections]
    sections_map = {s.get('className', ''): s for s in sections}

    # Node type from section classNames.
    if 'RandomSource' in sec_classes:
        node_spec['type'] = 'Source'
    elif 'JobSink' in sec_classes:
        node_spec['type'] = 'Sink'
    elif 'Storage' in sec_classes:
        # SPN Place: a Storage section holds tokens (capacity + drop rules).
        node_spec['type'] = 'Place'
    elif any(sc in ('Enabling', 'Timing', 'Firing') for sc in sec_classes):
        # SPN Transition: Enabling (input arcs/inhibitors) + Timing (modes,
        # firing distribution, servers) + Firing (output arcs).
        node_spec['type'] = 'Transition'
    elif 'Delay' in sec_classes:
        node_spec['type'] = 'Delay'
    elif 'ClassSwitch' in sec_classes:
        node_spec['type'] = 'ClassSwitch'
    elif 'Fork' in sec_classes:
        node_spec['type'] = 'Fork'
    elif 'Join' in sec_classes:
        node_spec['type'] = 'Join'
    elif any(sc in ('Server', 'PSServer') for sc in sec_classes):
        node_spec['type'] = 'Queue'

    # SPN nodes carry their token dynamics in dedicated sections; parse them
    # separately from the queueing sections below.
    if node_spec['type'] == 'Place':
        storage = sections_map.get('Storage')
        total_cap = None
        class_caps: Dict[str, int] = {}
        drop_rules: Dict[str, str] = {}
        if storage is not None:
            tc = _find_param(storage, 'totalCapacity')
            if tc is not None:
                v = tc.find('value')
                if v is not None and v.text:
                    try:
                        total_cap = int(v.text)
                    except ValueError:
                        pass
            cap_param = _find_param(storage, 'capacities')
            if cap_param is not None:
                for cls, sub in _walk_class_array(cap_param):
                    v = sub.find('value')
                    if v is not None and v.text:
                        try:
                            class_caps[cls] = int(v.text)
                        except ValueError:
                            pass
            drop_param = _find_param(storage, 'dropRules')
            if drop_param is not None:
                for cls, sub in _walk_class_array(drop_param):
                    v = sub.find('value')
                    if v is not None and v.text:
                        drop_rules[cls] = v.text.strip()
        node_spec['total_capacity'] = total_cap
        node_spec['class_capacities'] = class_caps
        node_spec['drop_rules'] = drop_rules
        return node_spec

    if node_spec['type'] == 'Transition':
        node_spec['modes'] = _parse_transition_modes(sections_map)
        return node_spec

    for section in sections:
        sc = section.get('className', '')

        if sc == 'RandomSource':
            node_spec['arrivals'] = _parse_service_array(section)

        elif sc in ('Server', 'PSServer', 'Delay'):
            services = _parse_service_array(section)
            if services:
                node_spec['services'] = services
            if sc == 'PSServer':
                node_spec['scheduling'] = 'PS'
            maxjobs = _find_param(section, 'maxJobs')
            if maxjobs is not None:
                value = maxjobs.find('value')
                if value is not None and value.text:
                    try:
                        nservers = int(value.text)
                        if nservers > 0:
                            node_spec['servers'] = nservers
                    except ValueError:
                        pass

        elif sc == 'Queue':
            size = _find_param(section, 'size')
            if size is not None:
                value = size.find('value')
                if value is not None and value.text:
                    try:
                        node_spec['capacity'] = int(value.text)
                    except ValueError:
                        pass
            # Scheduling get-strategy determines the queue discipline.
            for param in section.findall('parameter'):
                cp = param.get('classPath', '')
                if cp.endswith('PSStrategy'):
                    node_spec['scheduling'] = 'PS'
                elif cp.endswith('LCFSstrategy'):
                    node_spec['scheduling'] = 'LCFS'
                elif cp.endswith('FCFSstrategy'):
                    node_spec.setdefault('scheduling', 'FCFS')

        elif sc == 'Router':
            routing = _parse_routing_array(section)
            if routing:
                node_spec['routing'] = routing

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
        # JMT stores the closed population under the 'customers' attribute.
        try:
            pop = class_elem.get('customers', class_elem.get('population', '1'))
            class_spec['population'] = int(pop)
        except ValueError:
            class_spec['population'] = 1

    # Reference station/source. JMT names the attribute 'referenceSource' for
    # both open (the Source) and closed (the reference station) classes.
    ref_station = class_elem.get('referenceSource', class_elem.get('referenceStation', ''))
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

    # Max samples - support both 'samples' and 'max_samples' keys
    max_samples = ET.SubElement(params, 'maxSamples')
    max_samples.text = str(options.get('samples', options.get('max_samples', 1000000)))


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

    # Queue capacity (size parameter)
    # LINE uses Kendall notation where cap = K = total system capacity
    # JMT's "size" parameter represents total capacity K (-1 means infinite)
    capacity = -1  # Default: infinite capacity
    if hasattr(sn, 'cap') and sn.cap is not None and ist < len(sn.cap):
        cap_val = sn.cap[ist]
        if not np.isinf(cap_val):
            capacity = int(cap_val)

    size_param = ET.SubElement(queue_section, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = str(capacity)

    # Drop strategies - what to do when queue is full
    # For finite capacity queues (M/M/1/K), need to specify "drop" behavior
    drop_strat_param = ET.SubElement(queue_section, 'parameter')
    drop_strat_param.set('array', 'true')
    drop_strat_param.set('classPath', 'java.lang.String')
    drop_strat_param.set('name', 'dropStrategies')

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        ref_class = ET.SubElement(drop_strat_param, 'refClass')
        ref_class.text = class_name

        subparam = ET.SubElement(drop_strat_param, 'subParameter')
        subparam.set('classPath', 'java.lang.String')
        subparam.set('name', 'dropStrategy')

        value_elem = ET.SubElement(subparam, 'value')
        # Check if there's a drop rule defined in sn.droprule
        drop_text = 'drop'  # Default to drop for finite capacity
        if hasattr(sn, 'droprule') and sn.droprule is not None:
            if ist < sn.droprule.shape[0] and k < sn.droprule.shape[1]:
                drop_val = sn.droprule[ist, k]
                if drop_val == 0 or np.isnan(drop_val):
                    drop_text = 'drop'
                elif drop_val == 1:  # WAITQ
                    drop_text = 'waiting queue'
                else:
                    drop_text = 'drop'
        value_elem.text = drop_text

    # Service section (with number of servers)
    service_section = ET.SubElement(node, 'section')
    service_section.set('className', 'Server')

    # Number of servers (maxJobs parameter in Server section)
    # Use the maximum of nservers and lldscaling (for load-dependent models).
    # Load-dependent scaling with min(1:N,c) represents a c-server queue,
    # where max(lldscaling) = c.
    servers = 1
    if hasattr(sn, 'nservers') and ist < len(sn.nservers):
        srv_val = sn.nservers[ist]
        if np.isinf(srv_val):
            servers = -1  # Infinite servers
        else:
            servers = int(srv_val)

    # Check load-dependent scaling for effective number of servers
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        if ist < sn.lldscaling.shape[0]:
            effective_servers = int(np.max(sn.lldscaling[ist, :]))
            if servers > 0:  # Don't override infinite servers
                servers = max(servers, effective_servers)

    maxjobs_param = ET.SubElement(service_section, 'parameter')
    maxjobs_param.set('classPath', 'java.lang.Integer')
    maxjobs_param.set('name', 'maxJobs')
    value = ET.SubElement(maxjobs_param, 'value')
    value.text = str(servers)

    # Service strategies (in Server section)
    service_strat_param = ET.SubElement(service_section, 'parameter')
    service_strat_param.set('array', 'true')
    service_strat_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    service_strat_param.set('name', 'ServiceStrategy')

    for k in range(sn.nclasses):
        class_name = sn.classnames[k] if hasattr(sn, 'classnames') else f'Class{k}'

        # Get service rate
        rate = 0
        if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1]:
            rate = sn.rates[ist, k]
            if np.isnan(rate):
                continue

        mean = 1.0 / rate if rate > 0 else 0

        # Add refClass element
        ref_class = ET.SubElement(service_strat_param, 'refClass')
        ref_class.text = class_name

        subparam = ET.SubElement(service_strat_param, 'subParameter')
        subparam.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        subparam.set('name', 'ServiceTimeStrategy')

        dist_subparam = ET.SubElement(subparam, 'subParameter')
        dist_subparam.set('classPath', 'jmt.engine.random.Exponential')
        dist_subparam.set('name', 'Exponential')

        par_subparam = ET.SubElement(subparam, 'subParameter')
        par_subparam.set('classPath', 'jmt.engine.random.ExponentialPar')
        par_subparam.set('name', 'distrPar')

        lambda_param = ET.SubElement(par_subparam, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        value.text = str(rate) if rate > 0 else '1.0'


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


def _add_jsimg_regions(sim: ET.Element, model: Any, sn: Any) -> None:
    """
    Add finite capacity regions (blocking regions) to JSIMG.

    Creates blockingRegion elements for each FCR defined in the model.

    Args:
        sim: Parent sim XML element
        model: Network model (to access regions)
        sn: NetworkStruct

    References:
        MATLAB: matlab/src/solvers/JMT/@JMTIO/saveRegions.m
    """
    from ..sn.network_struct import DropStrategy as DSEnum

    # Get regions from model
    regions = []
    if hasattr(model, 'get_regions'):
        regions = model.get_regions()
    elif hasattr(model, 'regions'):
        regions = model.regions

    if not regions:
        return

    for r_idx, region in enumerate(regions):
        blocking_region = ET.SubElement(sim, 'blockingRegion')
        region_name = region.get_name() if hasattr(region, 'get_name') else f'FCRegion{r_idx + 1}'
        blocking_region.set('name', region_name)
        blocking_region.set('type', 'default')

        # 1. regionNode elements - nodes in this region
        region_nodes = region.nodes if hasattr(region, 'nodes') else []
        for node in region_nodes:
            node_name = node.get_name() if hasattr(node, 'get_name') else str(node)
            region_node = ET.SubElement(blocking_region, 'regionNode')
            region_node.set('nodeName', node_name)

        # 2. globalConstraint
        global_constraint = ET.SubElement(blocking_region, 'globalConstraint')
        global_max = region.global_max_jobs if hasattr(region, 'global_max_jobs') else -1
        global_constraint.set('maxJobs', str(global_max))

        # 3. globalMemoryConstraint
        global_mem_constraint = ET.SubElement(blocking_region, 'globalMemoryConstraint')
        global_max_mem = region.global_max_memory if hasattr(region, 'global_max_memory') else -1
        global_mem_constraint.set('maxMemory', str(global_max_mem))

        # Get classes from region or sn
        region_classes = region.classes if hasattr(region, 'classes') else []

        # 4. classConstraint elements
        for job_class in region_classes:
            class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
            class_max_jobs = region.get_class_max_jobs(job_class) if hasattr(region, 'get_class_max_jobs') else -1

            # Only write if not unbounded (-1)
            if class_max_jobs != -1:
                class_constraint = ET.SubElement(blocking_region, 'classConstraint')
                class_constraint.set('jobClass', class_name)
                class_constraint.set('maxJobsPerClass', str(class_max_jobs))

        # 5. classMemoryConstraint elements
        for job_class in region_classes:
            class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
            class_max_mem = region.get_class_max_memory(job_class) if hasattr(region, 'get_class_max_memory') else -1

            # Only write if not unbounded (-1)
            if class_max_mem != -1:
                class_mem_constraint = ET.SubElement(blocking_region, 'classMemoryConstraint')
                class_mem_constraint.set('jobClass', class_name)
                class_mem_constraint.set('maxMemoryPerClass', str(class_max_mem))

        # 6. dropRules elements - always write for each class
        for job_class in region_classes:
            class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
            drop_rule = region.get_drop_rule(job_class) if hasattr(region, 'get_drop_rule') else None

            drop_rules = ET.SubElement(blocking_region, 'dropRules')
            drop_rules.set('jobClass', class_name)

            # Determine if DROP or WAITQ
            if drop_rule is not None:
                # Check if it's DROP strategy
                is_drop = False
                if hasattr(drop_rule, 'name'):
                    is_drop = drop_rule.name == 'DROP'
                elif hasattr(drop_rule, 'value'):
                    is_drop = drop_rule.value == 1  # DROP = 1
                drop_rules.set('dropThisClass', 'true' if is_drop else 'false')
            else:
                drop_rules.set('dropThisClass', 'false')

        # 7. classSize elements (only if not default value of 1)
        for job_class in region_classes:
            class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
            class_size = region.get_class_size(job_class) if hasattr(region, 'get_class_size') else 1

            if class_size != 1:
                class_size_elem = ET.SubElement(blocking_region, 'classSize')
                class_size_elem.set('jobClass', class_name)
                class_size_elem.set('size', str(class_size))


__all__ = [
    'JSIMGInfo',
    'jmt2line',
    'jsim2line',
    'jmva2line',
    'qn2jsimg',
]
