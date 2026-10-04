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


def parse_jmt(filename: str) -> Dict[str, Any]:
    """
    Parse a JMT model file into a dictionary specification (no Network is built).

    Handles JSIMG (simulation) and JMVA (analytical) formats. Use jmt2line, or
    line_solver.io.JMT2LINE, for the Network that MATLAB's JMT2LINE returns.

    Args:
        filename: Path to JMT file (.jsimg or .jmva)

    Returns:
        Dictionary with network model specification

    References:
        MATLAB: matlab/src/io/JMT2LINE.m
    """
    ext = os.path.splitext(filename)[1].lower()

    if ext == '.jsimg':
        return parse_jsim(filename)
    elif ext == '.jmva':
        return parse_jmva(filename)
    else:
        # Try to detect format from file content
        try:
            tree = ET.parse(filename)
            root = tree.getroot()
            if root.tag == 'sim' or 'simType' in root.attrib:
                return parse_jsim(filename)
            elif root.tag == 'model' or 'modelType' in root.attrib:
                return parse_jmva(filename)
            else:
                line_warning('parse_jmt', f'Unknown JMT format in {filename}, trying JSIM')
                return parse_jsim(filename)
        except Exception as e:
            line_error('parse_jmt', f'Failed to parse JMT file: {e}')
            return {}


def parse_jsim(filename: str) -> Dict[str, Any]:
    """
    Parse a JSIM model file into a dictionary specification (no Network is built).

    JSIM files are XML-based simulation model definitions used by JMT. This is the
    parser line_solver.io.JSIM2LINE builds its Network from; use jsim2line, or
    JSIM2LINE, for the Network that MATLAB's JSIM2LINE returns.

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
        line_error('parse_jsim', f'XML parse error: {e}')
    except Exception as e:
        line_error('parse_jsim', f'Error importing JSIM: {e}')

    return result


def parse_jmva(filename: str) -> Dict[str, Any]:
    """
    Parse a JMVA model file into a coarse dictionary specification (no Network is built).

    JMVA files are XML-based analytical model definitions used by JMT. Use jmva2line,
    or line_solver.io.JMVA2LINE, for the Network that MATLAB's JMVA2LINE returns.

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
        line_error('parse_jmva', f'XML parse error: {e}')
    except Exception as e:
        line_error('parse_jmva', f'Error importing JMVA: {e}')

    return result


def jmt2line(filename: str, model_name: Optional[str] = None):
    """
    Import a JMT model file (.jmva, .jsim, .jsimg, .jsimw) as a LINE Network.

    Same semantics as MATLAB's JMT2LINE(filename, modelName); see
    line_solver.io.JMT2LINE. For the dictionary form use parse_jmt.

    References:
        MATLAB: matlab/src/io/JMT2LINE.m
    """
    from ...io import JMT2LINE
    return JMT2LINE(filename, model_name)


def jsim2line(filename: str, model_name: Optional[str] = None):
    """
    Import a JMT JSIM/JSIMG/JSIMW model file as a LINE Network.

    Same semantics as MATLAB's JSIM2LINE(filename, modelName); see
    line_solver.io.JSIM2LINE. For the dictionary form use parse_jsim.

    References:
        MATLAB: matlab/src/io/JSIM2LINE.m
    """
    from ...io import JSIM2LINE
    return JSIM2LINE(filename, model_name)


def jmva2line(filename: str, model_name: Optional[str] = None):
    """
    Import a JMT JMVA model file as a LINE Network.

    Same semantics as MATLAB's JMVA2LINE(filename, modelName); see
    line_solver.io.JMVA2LINE. For the dictionary form use parse_jmva.

    References:
        MATLAB: matlab/src/io/JMVA2LINE.m
    """
    from ...io import JMVA2LINE
    return JMVA2LINE(filename, model_name)


def qn2jsimg(model: Any, output_filename: Optional[str] = None,
             options: Optional[Dict] = None) -> str:
    """
    Export a Network model to JMT JSIMG format.

    Creates a JSIMG (JMT simulation) XML file from a Network model.

    Args:
        model: Network model or NetworkStruct
        output_filename: Optional output file path (default: temp file)
        options: SolverJMTOptions, a dict or an options object (samples, seed, ...); None for the defaults

    Returns:
        Path to the created JSIMG file

    References:
        MATLAB: matlab/src/io/QN2JSIMG.m
    """
    if output_filename is None:
        fd, output_filename = tempfile.mkstemp(suffix='.jsimg')
        os.close(fd)

    # As QN2JSIMG.m: the document is the one the JMT solver writes (JMTIO.writeJSIM), so what
    # this exports is what SolverJMT simulates and what JSIM2LINE reads back.
    from ..solvers.jmt.handler import _write_jsim_file, SolverJMTOptions
    if hasattr(model, 'getStruct'):
        sn, net = model.getStruct(), model
    else:
        sn, net = model, None
    _write_jsim_file(sn, output_filename, _jsimg_writer_options(options), model=net)
    return output_filename


def _jsimg_writer_options(options: Any):
    """SolverJMTOptions from None, a dict, or any options object carrying some of its fields."""
    from dataclasses import fields as _fields
    from ..solvers.jmt.handler import SolverJMTOptions
    if isinstance(options, SolverJMTOptions):
        return options
    names = [f.name for f in _fields(SolverJMTOptions)]
    if options is None:
        return SolverJMTOptions()
    if isinstance(options, dict):
        return SolverJMTOptions(**{k: v for k, v in options.items() if k in names})
    return SolverJMTOptions(**{k: getattr(options, k) for k in names if hasattr(options, k)})


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


def _param_int(section: ET.Element, name: str) -> Optional[int]:
    """Read a scalar integer <parameter name=...> of a section."""
    param = _find_param(section, name)
    if param is None:
        return None
    value = param.find('value')
    if value is None or not value.text:
        return None
    try:
        return int(float(value.text))
    except ValueError:
        return None


def _parse_class_switch_matrix(section: ET.Element) -> Dict[str, Dict[str, float]]:
    """Parse the ClassSwitch section into {from class: {to class: probability}}.

    The matrix is a per-class array of per-class rows, so the interleaved
    refClass/subParameter layout appears twice, once per index.
    """
    matrix: Dict[str, Dict[str, float]] = {}
    param = _find_param(section, 'matrix')
    if param is None:
        return matrix
    for from_class, row in _walk_class_array(param):
        entries: Dict[str, float] = {}
        for to_class, cell in _walk_class_array(row):
            value = cell.find('value')
            if value is not None and value.text:
                try:
                    entries[to_class] = float(value.text)
                except ValueError:
                    pass
        matrix[from_class] = entries
    return matrix


def _parse_join_required(section: ET.Element) -> Dict[str, int]:
    """Parse the Join section into {class: numRequired}, JMT's -1 meaning all."""
    required: Dict[str, int] = {}
    param = _find_param(section, 'JoinStrategy')
    if param is None:
        return required
    for class_name, strat in _walk_class_array(param):
        for sub in strat.findall('subParameter'):
            if sub.get('name') != 'numRequired':
                continue
            value = sub.find('value')
            if value is not None and value.text:
                try:
                    required[class_name] = int(float(value.text))
                except ValueError:
                    pass
    return required


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
    """Parse the per-class service/arrival distribution array of a
    RandomSource/Server/PSServer/Delay section into {className: distr_spec_or_None}.

    JMT names the distribution array 'ServiceStrategy' in RandomSource, Server and
    Delay sections but 'ServerStrategy' in a PSServer section, so both keys are
    accepted here; otherwise the per-class service rates of a processor-sharing
    queue are silently dropped (MATLAB/JAR read this array positionally and do
    not hit the mismatch)."""
    result = {}
    param = _find_param(section, 'ServiceStrategy')
    if param is None:
        param = _find_param(section, 'ServerStrategy')
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
        node_spec['csmatrix'] = _parse_class_switch_matrix(sections_map['ClassSwitch'])
    elif 'Fork' in sec_classes:
        node_spec['type'] = 'Fork'
        tasks = _param_int(sections_map['Fork'], 'jobsPerLink')
        if tasks is not None:
            node_spec['tasks_per_link'] = tasks
    elif 'Join' in sec_classes:
        node_spec['type'] = 'Join'
        node_spec['join_required'] = _parse_join_required(sections_map['Join'])
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


__all__ = [
    'JSIMGInfo',
    'jmt2line',
    'jsim2line',
    'jmva2line',
    'parse_jmt',
    'parse_jsim',
    'parse_jmva',
    'qn2jsimg',
]
